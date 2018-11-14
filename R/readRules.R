# library(data.table)
# library(plyr)
# library(dplyr)
# library(tidyr)
# library(readr)

#' Read network from BNET file
#' 
#' @param filename BNET file. The file is a tab delimited file with two columns,
#' `targets` and `factors`. `factors` contains the logic rule and `targets` the node
#' activated by the logic rule.
#'
#' @return CellNOpt network
#'
#' @author Luis Tobalina
readBNET <- function(filename){
	required_pcks = list("plyr","dplyr","tidyr","readr")
	
	if(!all(unlist(lapply(required_pcks,function(str) {require(str,character.only = TRUE)})))){
		print("the following packages need to be installed to use readBNET:")	
		print(unlist(required_pcks))
		print("Please, install the packages manually for this feature.")
		return(-1)
	}
			   
	warning("experimental BNET reader. Use with care (February 2018).")
	
	bnet <- data.table::fread(filename, header=T)
	
	# count number of `and` gates in each logic rule
	bnet <- bnet %>% rowwise() %>% 
		mutate(i_and_gates = (nchar(factors) - nchar(gsub("&", "", factors, perl=TRUE))))
	# count how many `and` gates have appeared up to the current rule
	bnet <- bnet %>% transform(i_and_gates = c(0, i_and_gates[-nrow(bnet)])) %>% 
		mutate(i_and_gates = cumsum(i_and_gates))
	# parse each logic rule
	sif <- bnet %>% rowwise() %>% 
		do(sif_df = build_sif_table_from_rule(.$factors, .$targets, last_and_num=.$i_and_gates)$sif_str) %>% 
		unnest() %>% 
		ungroup() 
	
	# use `readSIF()` to get the network ready for CellNOpt
	fh <- tempfile()
	write.table(sif$sif_df, file=fh,
				row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
	model <- readSIF(fh)
	
	return(model)
}

#' Read network from BND file
#' 
#' BND is a file format used by MaBoSS to store the boolean network definition.
#' The reader works if the logic for the activation of the node is stated withing the parameter `logic` of the node.
#' An example file can be found in `https://maboss.curie.fr/pub/example.bnd`.
#' 
#' @param filename BND file. 
#'
#' @return CellNOpt network
#' 
#' @author Luis Tobalina
readBND <- function(filename){
	required_pcks = list("plyr","dplyr","tidyr","readr")
	
	if(!all(unlist(lapply(required_pcks,function(str) {require(str,character.only = TRUE)})))){
		print("the following packages need to be installed to use readBND:")	
		print(unlist(required_pcks))
		print("Please, install the packages manually for this feature.")
		return(-1)
	}
	warning("experimental BND reader. Use with care (July 2018).")
	
	bnd <- read_file(filename)
	# remove all line breaks
	bnd <- gsub("\\n", "", bnd, perl=TRUE)
	# remove parenthesis enclosing a single node (which is composed of alphanumeric characters and "_")
	bnd <- gsub("\\((!?[[:alnum:]|_]*)\\)", "\\1", bnd, perl=TRUE)
	# get the information in bnet format 
	targets <- regmatches(bnd, gregexpr("(?<=Node ).*?(?= {)", bnd, perl = TRUE))[[1]]
	factors <- regmatches(bnd, gregexpr("(?<=logic = ).*?(?=;)", bnd, perl = TRUE))[[1]]
	bnet <- data.frame(targets, factors)
	
	# use `readBNET()` to get the network ready for CellNOpt
	fh <- tempfile()
	write.table(bnet, file=fh,
				row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
	model <- readBNET(fh)
}


# logic rule parser written following explanation of boolean expression evaluator in:
# https://unnikked.ga/how-to-build-a-boolean-expression-evaluator-518e9e068a65
# also of interest:
# https://www.strchr.com/expression_evaluator
# Unless the boolean expression is written in Disjunctive Normal Form (DNF), in order
# to get a correct `*.sif` file we will need to process the first interpretation we do of
# the tree into sif format.
#
#' Build a SIF table from a logic rule written in a string
#' 
#' @param rule_str String containing the rule to be parsed 
#' @param target Name of the node affected by the rule
#' @param last_and_num If the rule contains `and` gates, their numeration will start after the number provided here (default is 0)
#'
#' @return data.frame with the network structure derived from the rule.
#' The column `sif_str` contains the string that can be written to a file 
#' and then read with `readSIF()` in order to load a CellNOpt compatible network.
#' 
#' @examples
#' CellNOptR:::build_sif_table_from_rule("B & (C | D)", "A", last_and_num=2)
#' 
#' test_rule <- list()
#' test_rule[[1]] <- "AMP_ATP | (ATM & ATR) | HIF1 | !(EGFR | FGFR3)"
#' test_rule[[2]] <- "A & ((B | C) & !(D & E))"
#' test_rule[[3]] <- "A & B | C"
#' test_rule[[4]] <- "A & B & C"
#' test_rule[[5]] <- "A & (B | C)"
#' test_rule[[6]] <- "(A | B) & (C | D)"
#' test_rule[[7]] <- "!(C & D) | (E & F)"
#' test_rule[[8]] <- "(A | B) & (C | !D) & (E | F)"
#' parsed_rule <- list()
#' for (i in c(1:length(test_rule))){
#'   parsed_rule[[i]] <- CellNOptR:::build_sif_table_from_rule(test_rule[[i]], "T")
#' }
#' 
build_sif_table_from_rule <- function(rule_str, target, last_and_num=0) {
	required_pcks = list("plyr","dplyr","tidyr","readr")
	
	if(!all(unlist(lapply(required_pcks,function(str) {require(str,character.only = TRUE)})))){
		print("the following packages need to be installed to use readBND:")	
		print(unlist(required_pcks))
		print("Please, install the packages manually for this feature.")
		return(-1)
	}
	# remove all whitespaces
	rule_str <-  gsub(" ", "", rule_str)
	# split rule in its different atomic components
	tokens <- strsplit(rule_str, "(?=[\\|\\&\\(\\)!])", perl=T)[[1]]
	# the next set of variables will be used and/or modified by the inner functions
	current_token_pos <- 0
	symbol <- ""
	root <- ""
	negation_flag <- FALSE
	sif_str <- ""
	sif_list <- list()
	or_list <- list()
	and_list <- list()
	not_list <- list()
	sif_num <- 0
	or_num <- 0
	and_num <- last_and_num
	not_num <- 0
	
	# This function is used to retrieve the next token from the rule string
	get_next_token <- function(){
		current_token_pos <<- current_token_pos + 1
		#cat("\n", current_token_pos, tokens[current_token_pos], sep=" ")
		if (current_token_pos<=length(tokens)){
			symbol <<- tokens[current_token_pos]
			# apply de Morgan's Law if necessary
			if (negation_flag){
				if (symbol=="|"){
					symbol <<- "&"
				} else if (symbol=="&"){
					symbol <<- "|"
				}
			}
		} else {
			symbol <<- ""
		}
	}
	
	# helper function to parse each expression
	get_expression <- function(){
		get_term()
		while (symbol=="|") {
			left_part <- root
			get_term()
			right_part <- root
			
			or_num <<- or_num + 1
			or_list[[or_num]] <<- c(left_part, right_part)
			root <<- paste("or", or_num, sep="")
			
			sif_num <<- sif_num + 1
			sif_list[[sif_num]] <<- c(left_part, right_part, root)
		}
	}
	
	# helper function to parse each term
	get_term <- function(){
		get_factor()
		while (symbol=="&") {
			left_part <- root
			get_factor()
			right_part <- root
			
			and_num <<- and_num + 1
			and_list[[and_num]] <<- c(left_part, right_part)
			root <<- paste("and", and_num, sep="")
			
			sif_num <<- sif_num + 1
			sif_list[[sif_num]] <<- c(left_part, right_part, root)
		}
	}
	
	# helper function
	get_factor <- function(){
		get_next_token()
		if (symbol=="!") {
			# negate next expression
			negation_flag <<- !negation_flag
			get_factor()
			negation_flag <<- !negation_flag
			# revert changes to & and | symbols if they where acquired before the reset of the negation flag
			if (symbol=="|"){
				symbol <<- "&"
			} else if (symbol=="&"){
				symbol <<- "|"
			}
		} else if (symbol=="(") {
			get_expression()
			get_next_token()
		} else if (symbol==")") {
			# we don't care about ')'
		} else{
			if (negation_flag) {
				symbol <<- paste("!", symbol, sep="")
			}
			root <<- symbol
			get_next_token()
		}
	}
	
	# helper function to rename `and` or `or` gates (i.e. reset how they are numbered)
	rename_gates <- function(sif_df, type, last_num=0){
		if (!(type %in% c("and", "or"))){
			stop("type should be \"and\" or \"or\".")
		}
		type_gates <- unique(c(sif_df$node_in, sif_df$node_out))
		type_gates <- type_gates[grepl(paste("^", type, "\\d{1,}", sep=""), type_gates, perl=TRUE)]
		if (length(type_gates)>0){
			type_gates <- cbind(type_gates, paste("and", (c(1:length(type_gates))+last_num), sep=""))
			type_gates <- data.frame(type_gates)
			colnames(type_gates) <- c("old_name", "new_name")
			type_gates <- type_gates %>% mutate_if(is.factor, as.character)
			sif_df <- sif_df %>% 
				mutate(node_in = plyr::mapvalues(node_in, type_gates$old_name, type_gates$new_name, warn_missing=FALSE)) %>% 
				mutate(node_out = plyr::mapvalues(node_out, type_gates$old_name, type_gates$new_name, warn_missing=FALSE))
		}
		return(sif_df)
	}
	
	# helper function to simplify a chain of `and` or `or` gates
	simplify_gates <- function(sif_df, type){
		if (!(type %in% c("and", "or"))){
			stop("type should be \"and\" or \"or\".")
		}
		
		type_pattern <- paste("^", type, "\\d{1,}", sep="")
		while (any(grepl(type_pattern, sif_df$node_in, perl=TRUE) &
				   grepl(type_pattern, sif_df$node_out, perl=TRUE))){
			consecutive_type_gates <- sif_df %>% filter(grepl(type_pattern, sif_df$node_in, perl=TRUE), grepl(type_pattern, sif_df$node_out, perl=TRUE))
			sif_df <- sif_df %>% 
				mutate(node_in = plyr::mapvalues(node_in, consecutive_type_gates$node_in, consecutive_type_gates$node_out, warn_missing=FALSE)) %>% 
				mutate(node_out = plyr::mapvalues(node_out, consecutive_type_gates$node_in, consecutive_type_gates$node_out, warn_missing=FALSE)) %>% 
				filter(node_in!=node_out)
		}
		return(sif_df)
	}
	
	# helper function to interpret the parsed expression stored in the `sif_list` variable as a data.frame
	# `sif_list` contains the parsed expression in a tree like structure and this function takes care of making
	# several simplifications until a SIF format friendly version is reached (i.e. no concatenated `and` nodes and
	# no `or` nodes).
	interpret_sif_list <- function(sif_list){
		tree_df <- data.frame(matrix(unlist(sif_list), nrow=length(sif_list), byrow=T))
		colnames(tree_df) <- c("left_part", "right_part", "root")
		
		sif_df <- tree_df[,c("left_part", "root")]
		colnames(sif_df) <- c("node_in", "node_out")
		sif_df <- rbind(sif_df, setNames(tree_df[,c("right_part", "root")], names(sif_df)))
		# add target
		sif_df <- rbind(sif_df, data.frame(node_in=root, node_out=target))
		sif_df <- sif_df %>% mutate_if(is.factor, as.character)
		
		# simplify cascade of "and" and "or" gates
		sif_df <- simplify_gates(sif_df, "and")
		#sif_df <- simplify_gates(sif_df, "or")
		
		# we need to have the expression in disjunctive normal form (DNF)
		# find if there are or gates (node_in) connected to and gates (node_out)
		or_to_and_gates <- sif_df %>% filter(grepl("^or\\d{1,}", sif_df$node_in, perl=TRUE), grepl("^and\\d{1,}", sif_df$node_out, perl=TRUE))
		while (dim(or_to_and_gates)[1]>0){
			current_pair <- or_to_and_gates[1,]
			or_group_df <- sif_df %>% filter(grepl(current_pair$node_in, sif_df$node_out, perl=TRUE)) %>% 
				mutate(num_or = 1) %>% 
				mutate(num_or = cumsum(num_or)) %>% 
				mutate(new_node_out = paste(current_pair$node_out, node_out, num_or, sep="_")) %>% 
				mutate(old_node_out = node_out) %>% 
				mutate(node_out = new_node_out)
			and_group_df <- sif_df %>% filter(grepl(current_pair$node_out, sif_df$node_out, perl=TRUE) &
											  	!grepl(current_pair$node_in, sif_df$node_in, perl=TRUE)) %>% 
				mutate(new_node_out = paste(c(or_group_df$new_node_out), collapse=",")) %>% 
				mutate(new_node_out = strsplit(as.character(new_node_out), ",")) %>% 
				unnest(new_node_out) %>% 
				mutate(old_node_out = node_out) %>% 
				mutate(node_out = new_node_out)
			gate_group_df <- sif_df %>% filter(grepl(current_pair$node_out, sif_df$node_out, perl=TRUE) &
											   	grepl(current_pair$node_in, sif_df$node_in, perl=TRUE)) %>% 
				mutate(new_node_out = node_in) %>% 
				mutate(new_node_in = paste(c(or_group_df$new_node_out), collapse=",")) %>% 
				mutate(new_node_in = strsplit(as.character(new_node_in), ",")) %>% 
				unnest(new_node_in) %>% 
				mutate(old_node_in = node_in) %>% 
				mutate(old_node_out = node_out) %>% 
				mutate(node_in = new_node_in) %>% 
				mutate(node_out = new_node_out)
			root_group_df <- sif_df %>% filter(grepl(current_pair$node_out, sif_df$node_in, perl=TRUE)) %>% 
				mutate(new_node_in = current_pair$node_in) %>% 
				mutate(old_node_in = node_in) %>% 
				mutate(node_in = new_node_in)
			new_sif_df_part <- rbind(or_group_df[,c("node_in", "node_out")],
									 and_group_df[,c("node_in", "node_out")],
									 gate_group_df[,c("node_in", "node_out")],
									 root_group_df[,c("node_in", "node_out")])
			
			# if the or gate is connected to a different and gate than the one being processed now
			# we will rename it in those instances so that we don't lose these connections when filtering
			# the part of the dataframe that we are just about to change
			if (nrow(or_to_and_gates %>% filter(node_in==current_pair$node_in))>1) {
				duplicate_sif_df_part <- sif_df %>% filter((grepl(current_pair$node_in, sif_df$node_in, perl=TRUE) &
																!grepl(current_pair$node_out, sif_df$node_out, perl=TRUE)) |
																	(grepl(current_pair$node_in, sif_df$node_out, perl=TRUE)))
				or_num <<- or_num + 1
				new_or_gate_name <- paste("or", or_num, sep="")
				duplicate_sif_df_part <- duplicate_sif_df_part %>% 
					mutate(node_in = plyr::mapvalues(node_in, current_pair$node_in, new_or_gate_name, warn_missing=FALSE)) %>% 
					mutate(node_out = plyr::mapvalues(node_out, current_pair$node_in, new_or_gate_name, warn_missing=FALSE))
				sif_df <- sif_df %>% filter(!(grepl(current_pair$node_in, sif_df$node_in, perl=TRUE) &
											  	!grepl(current_pair$node_out, sif_df$node_out, perl=TRUE))) %>% 
					rbind(duplicate_sif_df_part)
			}
			
			# substitute the corresponding part of sif_df with the newly calculated df
			sif_df <- sif_df %>% filter(!((grepl(paste("^", current_pair$node_in, "$", sep=""), sif_df$node_out, perl=TRUE)) |
										  	(grepl(paste("^", current_pair$node_out, "$", sep=""), sif_df$node_out, perl=TRUE) &
										  	 	!grepl(paste("^", current_pair$node_in, "$", sep=""), sif_df$node_in, perl=TRUE)) |
										  	(grepl(paste("^", current_pair$node_out, "$", sep=""), sif_df$node_out, perl=TRUE) &
										  	 	grepl(paste("^", current_pair$node_in, "$", sep=""), sif_df$node_in, perl=TRUE)) |
										  	grepl(paste("^", current_pair$node_out, "$", sep=""), sif_df$node_in, perl=TRUE)))
			
			sif_df <- rbind(sif_df, new_sif_df_part)
			sif_df <- sif_df %>% distinct()
			
			# simplify cascade of "and" and "or" gates
			sif_df <- simplify_gates(sif_df, "and")
			sif_df <- simplify_gates(sif_df, "or")
			
			# update list of or to and gates
			or_to_and_gates <- sif_df %>% filter(grepl("^or\\d{1,}", sif_df$node_in, perl=TRUE), grepl("^and\\d{1,}", sif_df$node_out, perl=TRUE))
		}
		
		sif_df <- rename_gates(sif_df, type="and", last_num=last_and_num)
		
		# list or gates
		or_gates <- sif_df %>% filter(grepl("^or\\d{1,}", node_out, perl=TRUE)) %>% 
			group_by(node_out) %>% 
			summarise(or_members = paste(node_in, collapse=",")) %>% 
			mutate(root = node_out) %>% 
			select(root, or_members) %>% 
			ungroup()
		
		
		# substitute or gates by their inputs
		while (any(grepl("^or\\d{1,}", sif_df$node_in, perl=TRUE))) {
			sif_df <- sif_df %>% mutate(node_in = plyr::mapvalues(node_in, or_gates$root, or_gates$or_members, warn_missing=FALSE))
			sif_df <- sif_df %>% filter(!grepl("^or\\d{1,}", node_out, perl=TRUE))
			sif_df <- sif_df %>% mutate(node_in = strsplit(as.character(node_in), ",")) %>% unnest(node_in)
		}
		return(sif_df)
	}
	
	# helper function to write the column with the string for the SIF file
	write_sif <- function(sif_df){
		sif_df <- sif_df %>% mutate(sign1 = !grepl("^!", node_in, perl=TRUE)) %>% 
			mutate(sign2 = !grepl("^!", node_out, perl=TRUE)) %>% 
			mutate(sign = ifelse(sign1 & sign2, "1", "-1")) %>% 
			mutate(sif_str = paste(node_in, sign, node_out, sep="\t")) %>% 
			mutate(sif_str = gsub("^!", "", sif_str, perl=TRUE))
	}
	
	# After having defined all the helper functions, we parse the input expression.
	
	# if the rule only has one input (or a negated input), we output the sif format directly
	if (length(tokens)==1){
		sif_df <- data.frame(node_in=tokens, node_out=target)
		sif_df <- write_sif(sif_df)
		return(sif_df)
	} else if (length(tokens)==2){
		if (tokens[1]=="!"){
			sif_df <- data.frame(node_in=paste("!", tokens[2], sep=""), node_out=target)
			sif_df <- write_sif(sif_df)
			return(sif_df)
		}
	}
	
	# if there is more than one element, we parse the expression
	get_expression()
	sif_df <- interpret_sif_list(sif_list)
	sif_df <- write_sif(sif_df)
	
	#return(list(sif_df, sif_list))
	return(sif_df)
}
