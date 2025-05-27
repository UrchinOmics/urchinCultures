#functions and related helper code for x_0092 analysis


################################################################################
################################################################################
################################################################################



# Helper function: lookup LOC:GENENAME by common name
grab_feats <- function(seurat_obj, common_names) {
  all_features <- rownames(seurat_obj)
  matched <- sapply(common_names, function(name) {
    grep(paste0(":", name, "$"), all_features, value = TRUE)
  })
  unname(unlist(matched))
}


################################################################################
################################################################################
################################################################################



#' Find Feature Names Robustly and Optionally Log Actions
#'
#' Enhanced version of `grab_feats`. Finds feature names by trying exact match,
#' case-insensitive match, and matching by LOC ID prefix. Handles multiple matches
#' differently based on mode: interactively prompts user (with an 'All' option) or
#' uses specified non-interactive behavior ('first', 'all', 'none'). Can optionally
#' return a detailed log alongside the features.
#'
#' @param seurat_obj A Seurat object.
#' @param common_names A character vector of common gene names or LOC IDs to search for.
#' @param interactive Logical. If TRUE and multiple matches are found, prompt the user
#'                    to select the correct one (includes an "All" option). Defaults to FALSE.
#' @param multi_match_behavior Character. Specifies behavior when multiple matches are found
#'                             and `interactive = FALSE`. Options are "first" (default, take the
#'                             first match), "all" (take all unique matches), or "none" (take no matches).
#' @param return_log Logical. If TRUE, returns a list containing both the found features
#'                   and a detailed log tibble. If FALSE (default), returns only the character
#'                   vector of found features.
#'
#' @return If `return_log = FALSE` (default): A character vector of the unique feature names found.
#'         If `return_log = TRUE`: A list containing two elements:
#'         \describe{
#'           \item{`features`}{A character vector of the unique feature names found.}
#'           \item{`log`}{A tibble (data frame) detailing the process for each input name, including
#'                        query_name, status (e.g., "Single Match", "Multiple Matches (Interactive)"),
#'                        message (warnings issued), matches_found (potential matches), and
#'                        result_added (what was added to the `features` list).}
#'         }
#' @utility On a scale of 1-10: 9 (Very useful and robust gene lookup with flexible handling and optional logging)
#' @importFrom utils menu capture.output
#' @importFrom tibble tibble add_row
#' @importFrom dplyr bind_rows
#' @export
#' @examples
#' \dontrun{
#'   # Assuming 'seu' is a Seurat object
#'   # Default behavior: Get only the features vector
#'   features_only <- grab_feat_harder(seu, c("GeneA", "LOC12345", "AmbiguousGene"))
#'
#'   # Get both features and the log
#'   output_with_log <- grab_feat_harder(seu, c("GeneA", "NonExistent", "AmbiguousGene"),
#'                                       multi_match_behavior = "all", return_log = TRUE)
#'   found_features <- output_with_log$features
#'   log_df <- output_with_log$log
#'   print(log_df)
#'
#'   # --- Example: Managing the log externally ---
#'   # Initialize log outside function if needed
#'   # smelly_feat <- tibble::tibble(...)
#'
#'   # Run 1
#'   output1 <- grab_feat_harder(seu, c("GeneA", "GeneB"), return_log = TRUE)
#'   # if (!exists("smelly_feat")) smelly_feat <- output1$log else smelly_feat <- dplyr::bind_rows(smelly_feat, output1$log)
#'
#'   # Run 2
#'   output2 <- grab_feat_harder(seu, c("GeneC", "GeneD"), return_log = TRUE)
#'   # smelly_feat <- dplyr::bind_rows(smelly_feat, output2$log)
#'
#'   # Write log to file
#'   # if (exists("smelly_feat")) readr::write_tsv(smelly_feat, "feature_lookup_log.tsv")
#' }
grab_feat_harder <- function(seurat_obj,
                             common_names,
                             interactive = FALSE,
                             multi_match_behavior = c("first", "all", "none"),
                             return_log = FALSE) { # Added return_log argument
  
  # Ensure necessary packages are loaded if not using @importFrom
  # if (!requireNamespace("tibble", quietly = TRUE)) stop("Package 'tibble' needed.")
  # if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' needed.")
  
  # Validate multi_match_behavior argument
  multi_match_behavior <- match.arg(multi_match_behavior)
  
  all_features <- rownames(seurat_obj)
  results <- character()
  # Initialize log data frame (using tibble for better printing)
  log_collector <- list() # Collect rows as list elements first
  
  # Ensure common_names is not NULL or empty
  if (is.null(common_names) || length(common_names) == 0) {
    warning("Input 'common_names' is NULL or empty.")
    empty_log <- tibble::tibble(
      query_name = character(), status = character(), message = character(),
      matches_found = character(), result_added = character()
    )
    # Return based on return_log flag
    if (return_log) {
      return(list(features = character(0), log = empty_log))
    } else {
      return(character(0))
    }
  }
  
  for (name in common_names) {
    # Initialize log entry parts for this name
    current_query <- name
    current_status <- "Pending"
    current_message <- NA_character_
    current_matches_str <- NA_character_
    current_result_str <- NA_character_
    
    if (is.null(name) || name == "") {
      current_status <- "Skipped"
      current_message <- "Input name was NULL or empty."
      warning(current_message)
      # Add entry to log collector before skipping
      log_collector[[length(log_collector) + 1]] <- tibble::tibble(
        query_name = current_query, status = current_status, message = current_message,
        matches_found = current_matches_str, result_added = current_result_str
      )
      next
    }
    # Basic escaping for regex special characters in the input name
    safe_name <- gsub("([.^$*+?()[{|\\-])", "\\\\\\1", name)
    
    # --- Try matching ---
    matches <- character(0) # Ensure matches is defined
    match_step <- "Initial"
    
    # Try exact match first (suffix)
    matches <- grep(paste0(":", safe_name, "$"), all_features, value = TRUE)
    if (length(matches) > 0) match_step <- "Exact Suffix"
    
    # Relax to case-insensitive match if needed (suffix)
    if (length(matches) == 0) {
      matches <- grep(paste0(":", safe_name, "$"), all_features, value = TRUE, ignore.case = TRUE)
      if (length(matches) > 0) match_step <- "Case-insensitive Suffix"
    }
    
    # Try LOC ID prefix match if still nothing
    if (length(matches) == 0) {
      # Check if name looks like a LOC ID before trying prefix match
      if (grepl("^LOC[0-9]+$", name, ignore.case = TRUE)) {
        matches <- grep(paste0("^", safe_name, ":"), all_features, value = TRUE)
        if (length(matches) > 0) match_step <- "LOC ID Prefix"
      }
    }
    # --- End Matching ---
    
    
    # Handle results
    if (length(matches) == 0) {
      current_status <- "No Match"
      current_message <- paste("No match found for:", name)
      warning(current_message)
      
    } else if (length(matches) == 1) {
      current_status <- paste0("Single Match (", match_step, ")")
      current_matches_str <- matches[1]
      current_result_str <- matches[1]
      results <- c(results, matches[1])
      
    } else { # Multiple matches
      unique_matches <- unique(matches)
      current_matches_str <- paste(unique_matches, collapse=", ")
      multi_match_msg <- paste("Multiple unique matches found:", current_matches_str)
      # Add context to warning
      warning(paste("For query '", name, "': ", multi_match_msg, sep=""))
      
      if (interactive && base::interactive()) { # Handle interactive mode
        current_status <- "Multiple Matches (Interactive)"
        cat("\nMultiple unique matches found for:", name, "\n")
        # Add "All of the above" option to the menu choices
        menu_choices <- c(unique_matches, "All of the above")
        selected_index <- utils::menu(menu_choices, title = "Select the correct feature(s):")
        
        if (selected_index == 0) { # User cancelled
          current_message <- paste(multi_match_msg, "-> No selection made by user.", sep=" ")
          warning(paste("No selection made for:", name))
          current_result_str <- NA_character_ # Indicate nothing added
        } else if (selected_index == length(menu_choices)) { # User selected "All"
          current_message <- paste(multi_match_msg, "-> User selected 'All'.", sep=" ")
          results <- c(results, unique_matches)
          current_result_str <- paste(unique_matches, collapse=", ") # Log all added
          message("Selected all matches for: ", name)
        } else { # User selected a specific match
          selected_match <- unique_matches[selected_index]
          current_message <- paste(multi_match_msg, "-> User selected '", selected_match, "'.", sep="")
          results <- c(results, selected_match)
          current_result_str <- selected_match # Log the single one added
        }
        
      } else { # Handle non-interactive mode based on multi_match_behavior
        current_status <- paste0("Multiple Matches (Non-interactive: ", multi_match_behavior, ")")
        current_message <- paste(multi_match_msg, "-> Applying '", multi_match_behavior, "' behavior.", sep="")
        message(current_message) # Inform user via message
        
        switch(multi_match_behavior,
               first = {
                 results <- c(results, unique_matches[1])
                 current_result_str <- unique_matches[1]
               },
               all = {
                 results <- c(results, unique_matches)
                 current_result_str <- paste(unique_matches, collapse=", ")
               },
               none = {
                 current_result_str <- NA_character_ # Indicate nothing added
                 # Do nothing to results
               }
        )
      }
    } # End handling multiple matches
    
    # Add entry to log collector
    log_collector[[length(log_collector) + 1]] <- tibble::tibble(
      query_name = current_query, status = current_status, message = current_message,
      matches_found = current_matches_str, result_added = current_result_str
    )
    
  } # End loop over names
  
  # Combine log rows into a final tibble
  test_log_df <- dplyr::bind_rows(log_collector)
  
  # --- Return based on return_log flag ---
  if (return_log) {
    # Return list containing unique features and the log
    return(list(
      features = unique(results),
      log = test_log_df
    ))
  } else {
    # Default: return only the unique feature vector
    return(unique(results))
  }
}




################################################################################
################################################################################
################################################################################



match_feat_by_LOC <- function(obj, loc_ids) {
  all_feats <- rownames(obj)
  matched_feats <- all_feats[grepl(
    pattern = paste0("^(", paste(loc_ids, collapse = "|"), ")(:|$)"),
    x = all_feats
  )]
  return(matched_feats)
}




################################################################################
################################################################################
################################################################################


annotate_percent_ribo <- function(seurat_obj, ribo_genes) {
  message("Annotating object: ", substitute(seurat_obj))
  
  # make sure RNA assay is active and layer exists
  DefaultAssay(seurat_obj) <- "RNA"
  available_layers <- Layers(seurat_obj[["RNA"]])
  
  if (!"counts" %in% available_layers) {
    stop("No 'counts' layer found in RNA assay.")
  }
  
  ribo_genes <- intersect(ribo_genes, rownames(seurat_obj))
  if (length(ribo_genes) == 0) {
    warning("No ribosomal genes matched rownames in object.")
    return(seurat_obj)
  }
  
  counts <- GetAssayData(seurat_obj, layer = "counts")
  ribo_counts <- Matrix::colSums(counts[ribo_genes, , drop = FALSE])
  total_counts <- Matrix::colSums(counts)
  seurat_obj[["percent.ribo"]] <- ribo_counts / total_counts * 100
  
  return(seurat_obj)
}




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################






# --- Object Inspection ---

#' Describe Seurat Object Summary (v5+ focused with single verbose option)
#'
#' Prints a detailed summary of a Seurat object's structure, contents,
#' and processing history. A single verbose option controls the level of detail,
#' providing as much non-data matrix information as possible when TRUE.
#' Focused on Seurat v5+ object structure.
#'
#' @param obj A Seurat object.
#' @param name The name of the object variable (automatically captured).
#' @param n_meta Integer, number of metadata columns to preview in non-verbose mode (default 10).
#' @param n_levels Integer, number of factor levels to preview in non-verbose mode (default 5).
#' @param verbose Logical, if TRUE, prints verbose details for commands, reductions, graphs, and metadata, including sampled cell/feature names (default FALSE).
#'
#' @return Invisibly returns NULL. Output is printed to the console.
#' @utility On a scale of 1-10: 9 (Very comprehensive for understanding modern object state and history for analysis, controlled by verbose)
#' @importFrom Seurat Assays DefaultAssay Reductions Embeddings VariableFeatures Idents
#' @importFrom SeuratObject Layers # Relying directly on Layers from SeuratObject
#' @importFrom utils head tail packageVersion
#' @export
#' @examples
#' \dontrun{
#'   # Assuming 'seu' is a Seurat object (preferably v5+) and SeuratObject is loaded
#'   describe_seurat(seu) # Default concise summary
#'   describe_seurat(seu, verbose = TRUE) # Print maximum non-data matrix details
#' }
describe_seurat <- function(obj, name = deparse(substitute(obj)),
                            n_meta = 10, n_levels = 5,
                            verbose = FALSE) {
  
  cat("========================================\n")
  cat("====== Summary for Seurat Object:", name, "======\n")
  cat("========================================\n")
  cat("Class:", class(obj)[1], "\n")
  cat("Number of Cells:", ncol(obj), "\n")
  cat("Number of Features (Total):", nrow(obj), "\n")
  cat("----------------------------------------\n")
  
  # Software Versions (Helpful for debugging/documentation lookup)
  # Convert package_version object to character for cat()
  seurat_ver <- tryCatch(as.character(utils::packageVersion("Seurat")), error=function(e) "Unknown")
  seurat_object_ver <- tryCatch(as.character(utils::packageVersion("SeuratObject")), error=function(e) "Unknown")
  cat("Software Versions: Seurat v", seurat_ver, ", SeuratObject v", seurat_object_ver, "\n", sep="")
  cat("----------------------------------------\n")
  
  
  # Basic Identity Information
  cat("Active Assay:", Seurat::DefaultAssay(obj), "\n")
  active_idents <- tryCatch({ unique(Seurat::Idents(obj)) }, error=function(e) "[Error getting idents]")
  if (!identical(active_idents, "[Error getting idents]")) {
    cat("Active Ident (first 5):", paste(utils::head(active_idents, n = 5), collapse = ", "), if(length(active_idents) > 5) "..." else "", "\n")
    cat("Total Active Ident Levels:", length(active_idents), "\n")
  } else {
    cat("Active Ident: [Error getting idents]\n")
  }
  cat("----------------------------------------\n")
  
  # Sampled Names (Verbose Only)
  if (verbose) {
    cat("Sampled Cell Names (100 random):", "\n")
    # Check if ncol > 0 before sampling
    if (ncol(obj) > 0) {
      sampled_cells <- sample(colnames(obj), min(100, ncol(obj)))
      print(sampled_cells)
    } else {
      cat("  Object has 0 cells.\n")
    }
    cat("----------------------------------------\n")
    
    cat("Sampled Feature Names (100 random):", "\n")
    # Check if nrow > 0 before sampling
    if (nrow(obj) > 0) {
      sampled_features <- sample(rownames(obj), min(100, nrow(obj)))
      print(sampled_features)
    } else {
      cat("  Object has 0 features.\n")
    }
    cat("----------------------------------------\n")
  }
  
  
  # Processing History (@commands)
  commands_list <- obj@commands # Get the full list of commands
  commands_run <- names(commands_list)
  cat("Processing History (@commands):", length(commands_run), "commands recorded.\n")
  if (verbose) {
    if (length(commands_run) > 0) {
      # Print the full list of commands (non-data matrix parameters)
      print(commands_list)
    } else {
      cat("  No commands recorded in @commands slot.\n")
    }
  } else { # Concise summary
    if (length(commands_run) > 0) {
      cat("  Command Names (showing first 5 and last 5 if >10):\n")
      if (length(commands_run) > 10) {
        print(utils::head(commands_run, n = 5))
        print(utils::tail(commands_run, n = 5))
      } else {
        print(commands_run)
      }
      cat("  Inspect full command parameters with `", name, "@commands$CommandName`\n", sep="")
      cat("  Parameters for resulting reductions/graphs are stored within `[['ReductionName']]`/`[['GraphName']]` objects.\n")
    } else {
      cat("  No commands recorded in @commands slot.\n")
    }
  }
  cat("----------------------------------------\n")
  
  # Integration Status (Check command history for common integration steps)
  integration_commands <- c("IntegrateData", "RunHarmony", "RunFastMNN", "FindIntegrationAnchors", "RunAzimuth") # Added RunAzimuth
  used_integration_cmd <- intersect(commands_run, integration_commands)
  if (length(used_integration_cmd) > 0) {
    cat("Integration/Reference Mapping Status: Likely integrated or mapped (found commands: ", paste(used_integration_cmd, collapse = ", "), ")\n", sep="")
  } else {
    cat("Integration/Reference Mapping Status: Unlikely integrated/mapped via standard Seurat/common methods.\n")
  }
  cat("----------------------------------------\n")
  
  # QC Filtering Summary (Check command history for Subset)
  subset_commands <- commands_list[grep("^Subset", commands_run)]
  if (length(subset_commands) > 0) {
    cat("QC Filtering Status: Subset command(s) found in history.\n")
    for (sub_cmd_name in names(subset_commands)) {
      sub_cmd <- subset_commands[[sub_cmd_name]]
      filter_string <- tryCatch(sub_cmd$subset, error = function(e) NULL) # $subset stores the string filter
      filter_name_val <- tryCatch({ # Subset can also use subset.name/subset.value
        list(name = sub_cmd$subset.name, value = sub_cmd$subset.value)
      }, error = function(e) NULL)
      
      cat("  - Command '", sub_cmd_name, "': ", sep="")
      if (!is.null(filter_string)) {
        cat("Parameters: `", filter_string, "`\n", sep="")
      } else if (!is.null(filter_name_val) && !is.null(filter_name_val$name) && !is.null(filter_name_val$value)) {
        cat("Parameters: `", filter_name_val$name, " ", filter_name_val$value, "`\n", sep="")
      } else {
        cat("Parameters: Details unavailable\n")
      }
    }
    cat("  Review command history for exact parameters used.\n")
  } else {
    cat("QC Filtering Status: No standard Subset command found in history.\n")
  }
  cat("----------------------------------------\n")
  
  
  # Assay information
  assays_present <- Seurat::Assays(obj)
  cat("Assays Present:", paste(assays_present, collapse = ", "), "\n")
  
  for (assay in assays_present) {
    cat("  --- Assay:", assay, "---\n")
    assay_obj <- obj[[assay]]
    cat("    Features:", nrow(assay_obj), "\n")
    
    # Simplified Layer handling: Rely on SeuratObject::Layers
    cat("    Layers:")
    layers_present <- tryCatch({
      SeuratObject::Layers(assay_obj)
    }, error = function(e) {
      return(NULL) # Return NULL if getting layers fails
    })
    
    if (!is.null(layers_present) && length(layers_present) > 0) {
      # Print the names if successful and not empty
      cat(" ", paste(layers_present, collapse = ", "), "\n")
    } else {
      # Indicate if Layers() failed or returned empty/NULL
      cat(" None or Unknown\n")
    }
    
    # Variable features (Verbose Only)
    if (verbose) {
      var_feat <- tryCatch(Seurat::VariableFeatures(obj, assay = assay), error = function(e) NULL)
      if (!is.null(var_feat)) {
        cat("    Variable Features Count:", length(var_feat), "\n")
        # Optional: Could sample variable feature names here if needed
      } else {
        cat("    Variable Features Count: 0 or not found.\n")
      }
    }
  }
  cat("----------------------------------------\n")
  
  # Reduction information
  reductions_present <- Seurat::Reductions(obj)
  cat("Reductions Present:", paste(reductions_present, collapse = ", "), "\n")
  for (reduc in reductions_present) {
    cat("  --- Reduction:", reduc, "---\n")
    reduc_obj <- obj[[reduc]]
    # Embeddings access might still be slow for very large matrices, but ncol is fast
    embeddings_dim <- tryCatch({
      ncol(Seurat::Embeddings(reduc_obj))
    }, error = function(e) {
      return("[Error]")
    })
    cat("    Dimensions:", embeddings_dim, "\n")
    
    if (verbose) {
      cat("    Type:", class(reduc_obj)[1], "\n")
      # Try to get variance explained for PCA
      if (class(reduc_obj)[1] == "PCA") {
        var_exp <- tryCatch({ (reduc_obj@stdev)^2 }, error = function(e) NULL)
        if (!is.null(var_exp)) {
          cat("    Variance Explained (first 10):", paste(head(var_exp, 10), collapse = ", "), if(length(var_exp)>10) "..." else "", "\n")
        }
      }
      # Note about parameters being in commands or graph
      cat("    Full parameters often in @commands or related Graph object's @params slot.\n")
      # Optional: Could try to print reduc_obj@misc if it exists, but it's not standardized
      # misc_params <- tryCatch(reduc_obj@misc, error = function(e) NULL)
      # if (!is.null(misc_params)) { cat("    @misc slot:\n"); print(misc_params) }
    }
  }
  cat("----------------------------------------\n")
  
  # Graph information
  graphs_present <- names(obj@graphs)
  cat("Graphs Present:", paste(graphs_present, collapse = ", "), "\n")
  if (verbose) {
    for (graph_name in graphs_present) {
      cat("  --- Graph:", graph_name, "---\n")
      graph_obj <- obj[[graph_name]]
      cat("    Type:", class(graph_obj)[1], "\n") # SNN, KNN etc.
      
      # Try to get input details
      assay_used <- tryCatch(graph_obj@assay.used, error=function(e) "Unknown")
      cat("    Assay Used:", assay_used, "\n")
      # Access reduction/dims from @dim.reducs slot
      dims_details <- tryCatch(graph_obj@dim.reducs, error=function(e) list()) # Should be a named list
      if (length(dims_details) > 0) {
        # Iterate through entries in @dim.reducs (usually just one or two)
        for (reduc_name in names(dims_details)) {
          input_dims <- dims_details[[reduc_name]]
          cat("    Input Reduction & Dims:", reduc_name, "(Dims:", paste(input_dims, collapse=", "), ")\n")
        }
      } else {
        cat("    Input Reduction & Dims: Unknown or not recorded in @dim.reducs\n")
      }
      
      # Try to get parameters
      params <- tryCatch(graph_obj@params, error=function(e) NULL)
      if (!is.null(params) && length(params) > 0) {
        cat("    Parameters (@params slot):\n")
        print(params) # Print list of parameters
      } else {
        cat("    Parameters: None recorded in @params slot.\n")
      }
    }
  }
  cat("----------------------------------------\n")
  
  
  # Metadata information
  meta_cols <- colnames(obj[[]])
  cat("Metadata Columns:", length(meta_cols), "total.\n")
  
  # Suggest likely columns for biological context (Always shown)
  likely_bio_cols <- intersect(c("orig.ident", "group", "condition", "sample", "batch", "samplesource", "time", "timepoint"), tolower(meta_cols))
  if (length(likely_bio_cols) > 0) {
    cat("  Likely Biological Context Columns (based on names):", paste(meta_cols[match(likely_bio_cols, tolower(meta_cols))], collapse = ", "), "\n")
  } else {
    cat("  No common biological context columns found based on names (e.g., orig.ident, group, batch).\n")
  }
  
  
  if (verbose) {
    cat("  Column Names:\n")
    print(meta_cols) # Print all names
    
    cat("  Summary for Each Column:\n")
    for (col_name in meta_cols) {
      cat("  --- Metadata:", col_name, "---\n")
      # Use drop=TRUE to get a vector, which is better for summary functions
      col_data <- tryCatch(obj[[col_name, drop = TRUE]], error = function(e) "[Error accessing data]")
      
      if (identical(col_data, "[Error accessing data]")) {
        cat("    Class: [Error accessing data]\n")
        cat("    Summary: [Error accessing data]\n")
        next # Skip to next column if data access failed
      }
      
      
      cat("    Class:", class(col_data)[1], "\n")
      na_count <- sum(is.na(col_data))
      if (na_count > 0) {
        cat("    NA count:", na_count, "\n")
      }
      
      cat("    Summary:\n")
      # Use different summaries based on class
      if (is.numeric(col_data) || is.integer(col_data) || is.logical(col_data)) {
        print(summary(col_data))
      } else if (is.factor(col_data) || is.character(col_data)) {
        val_counts <- table(col_data, useNA = "ifany") # Include NA in table
        if (length(val_counts) > 0) {
          # Limit table output for very high cardinality
          max_levels_print = 50 # Define a limit for printing levels
          cat("    Levels/Values (", length(val_counts), " total):\n", sep="")
          if (length(val_counts) > max_levels_print) {
            print(head(val_counts, max_levels_print))
            cat("      ... and", length(val_counts) - max_levels_print, "more values.\n")
          } else {
            print(val_counts) # Print full table
          }
        } else {
          cat("    Levels/Values: None\n")
        }
      } else {
        # For other types, just print head
        cat("    Head of data:\n")
        print(utils::head(col_data))
      }
    }
    
  } else { # Non-verbose metadata output (concise)
    cat("  Column Names (showing first ", min(n_meta, length(meta_cols)), "):\n", sep="")
    print(utils::head(meta_cols, n = n_meta))
    
    # Preview levels of key factors
    key_factors <- c("orig.ident", grep("_clusters", meta_cols, value = TRUE), grep("_snn_res", meta_cols, value = TRUE))
    key_factors <- intersect(key_factors, meta_cols) # Keep only existing ones
    if (length(key_factors) > 0) {
      cat("  Preview of Key Metadata Factor Levels (showing first ", n_levels, "):\n", sep="")
      for (factor_col in key_factors) {
        # Use tryCatch here too, accessing columns for many cells can sometimes be slow
        levels_present <- tryCatch({
          unique(obj[[factor_col, drop=TRUE]])
        }, error = function(e) {
          return("[Error getting levels]")
        })
        
        if (!identical(levels_present, "[Error getting levels]")) {
          cat("  - ", factor_col, " (", length(levels_present), " levels): ",
              paste(utils::head(levels_present, n = n_levels), collapse=", "),
              if(length(levels_present) > n_levels) "..." else "", "\n", sep="")
        } else {
          cat("  - ", factor_col, ": [Error getting levels]\n", sep="")
        }
      }
    }
  }
  cat("========================================\n\n")
  invisible(NULL)
}





################################################################################
################################################################################
################################################################################


# Helper to sanitize filenames based on gene/tx annotations
cleanup <- function(name) {
  gsub("[/:\\\\\\s]", "_", name)  # Replace /, :, \, and whitespace with underscores
}


################################################################################
################################################################################
################################################################################
#need to test

get_limits <- function(seurat_obj, reduction = "umap") {
  available <- names(seurat_obj@reductions)
  
  if (is.null(reduction)) {
    stop("Please specify a reduction. Available options: ", paste(available, collapse = ", "))
  }
  
  reduction <- match.arg(reduction, choices = available)
  
  coords <- Embeddings(seurat_obj, reduction)
  
  list(
    xlims = range(coords[, 1], na.rm = TRUE),
    ylims = range(coords[, 2], na.rm = TRUE)
  )
}



################################################################################
################################################################################
################################################################################

#paired functions to find and set axis limits for plotting seurat objects. 
#call limited() on a seurat object and optionally specify a reduction. (default is umap)
#call limitless() with a plotting function to set axis limits which are appropriate for the data
#examples
# limitless(DimPlot(limited(obj), group.by = "orig.ident"))
# DimPlot(limited(obj),  group.by = "orig.ident") + limitless ()
# limited(obj)
# # DimPlot(obj) + limitless()
# # <any other plot> + limitless()

# Private environment to store axis limits
.limit_env <- new.env(parent = emptyenv())

# limited(): store x/y axis limits from any 2D reduction
limited <- function(seurat_obj, reduction = "umap") {
  available <- names(seurat_obj@reductions)
  if (!reduction %in% available) {
    stop("Reduction '", reduction, "' not found in object. Available reductions: ", paste(available, collapse = ", "))
  }
  coords <- Embeddings(seurat_obj, reduction)
  .limit_env$xlims <- range(coords[, 1], na.rm = TRUE)
  .limit_env$ylims <- range(coords[, 2], na.rm = TRUE)
  return(seurat_obj)
}

# limitless(): apply the stored limits to a ggplot-like plot
limitless <- function(p) {
  if (is.null(.limit_env$xlims) || is.null(.limit_env$ylims)) {
    stop("No limits set. Call limited(seurat_obj, reduction) first.")
  }
  p + xlim(.limit_env$xlims) + ylim(.limit_env$ylims)
}

################################################################################
################################################################################
################################################################################
#' MultiDimPlot
#'
#' Systematically generate and save `DimPlot` visualizations from a Seurat object
#' across multiple metadata fields and 2D dimension combinations.
#'
#' @param seu A single Seurat object.
#' @param reduction The name of the dimensional reduction to use (e.g., "umap", "pca").
#' @param ndims Number of dimensions to consider when building 2D plots (default = 2).
#' @param metadata_fields A character vector of metadata fields (columns in `seu@meta.data`) to use for coloring the plots.
#' @param output_dir Base directory where plots will be saved.
#' @param metadata_subdirectory Logical. If TRUE, saves each metadata field into its own subdirectory (default = TRUE).
#' @param return_plots Logical. If TRUE, returns a named list of generated ggplot objects (default = TRUE).
#' @param ... Additional arguments forwarded to `DimPlot()`, such as `label = TRUE`.
#'
#' @details
#' For a given Seurat object and reduction:
#' \itemize{
#'   \item All 2D combinations of the first `ndims` dimensions are plotted.
#'   \item For each metadata field specified:
#'     \itemize{
#'       \item A directory is created (optionally by metadata field) inside `output_dir`.
#'       \item A `DimPlot` is created colored by that metadata field.
#'       \item The plot title includes the Seurat object name, reduction, metadata field, and dimensions plotted.
#'       \item The plot is saved as a PNG file sized 14x7 inches.
#'       \item (Optional) The plot is returned inside a list.
#'     }
#' }
#'
#' The function automatically generates informative filenames and titles for all plots.
#'
#' @return A list of ggplot objects (if `return_plots = TRUE`), otherwise returns `NULL` invisibly.
#'
#' @import ggplot2
#' @import Seurat
#' @export
#'
#' @examples
#' MultiDimPlot(
#'   seu = my_seurat_object,
#'   reduction = "umap",
#'   ndims = 3,
#'   metadata_fields = c("orig.ident", "predicted.celltype"),
#'   output_dir = "plots/DimPlots",
#'   metadata_subdirectory = TRUE
#' )

MultiDimPlot <- function(seu,
                         reduction,
                         ndims,
                         metadata_fields,
                         output_dir,
                         metadata_subdirectory = TRUE,
                         return_plots = TRUE,
                         ...) {
  # capture the name of the Seurat object as a string
  obj_name <- deparse(substitute(seu))
  
  # 1. ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 2. build all 2-dimensional combos
  dims_combos <- combn(seq_len(ndims), 2, simplify = FALSE)
  
  # prepare list to store plots (only if needed)
  if (return_plots) {
    plots <- list()
  }
  
  # 3. loop over each metadata field
  for (field in metadata_fields) {
    # decide subdirectory for this field
    sub_dir <- if (metadata_subdirectory) {
      d <- file.path(output_dir, field)
      if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
      d
    } else {
      output_dir
    }
    
    # 4. for each dimension pair, make & save the plot
    for (dims in dims_combos) {
      dim_str    <- paste(dims, collapse = "_")
      
      # now include object name in title
      plot_title <- paste0(obj_name, " • ", reduction, 
                           " • ", field, 
                           " • dims ", dims[1], " vs ", dims[2])
      
      # and in the filename
      file_name  <- paste0("MultiDimPlot_",
                           obj_name, "_",
                           field, "_dims", dim_str, ".png")
      file_path  <- file.path(sub_dir, file_name)
      
      # create the ggplot object, forwarding any extra args
      p <- DimPlot(seu,
                   reduction = reduction,
                   group.by  = field,
                   dims      = dims,
                   ...) +
        ggtitle(plot_title) +
        theme(aspect.ratio = 1)
      
      # save it at 14×7 inches
      ggsave(filename = file_path,
             plot     = p,
             width    = 14,
             height   = 7,
             units    = "in")
      
      # store it if requested
      if (return_plots) {
        plots[[ paste(obj_name, field, dim_str, sep = "_") ]] <- p
      }
    }
  }
  
  # 5. return (or not)
  if (return_plots) {
    return(plots)
  } else {
    invisible(NULL)
  }
}

################################################################################
################################################################################
################################################################################



#good goi function (split umap expression data in ref-query projection space)
#good_goi, unlike dense goi, does not do any smoothing or computation of expression density
good_goi <- function(goi, ref, query,
                     ref_reduction,
                     query_reduction,
                     filename_prefix = "GoodGoiPlot_",
                     layout = "horizontal",
                     xlim = NULL,
                     ylim = NULL,
                     grey_zero = TRUE,
                     viridis_option = "D") {
  
  layout_op <- switch(layout,
                      vertical = `/`,
                      horizontal = `+`,
                      stop("layout must be 'horizontal' or 'vertical'"))
  
  get_plot_df <- function(obj, reduction_name, dimnames) {
    embed <- Embeddings(obj, reduction = reduction_name)
    exprs <- FetchData(obj, vars = goi)[, 1]
    df <- data.frame(x = embed[, 1], y = embed[, 2], expr = exprs)
    df <- df[order(df$expr), ]
    names(df)[1:2] <- dimnames
    return(df)
  }
  
  ref_dims <- paste0(ref_reduction, 1:2)
  query_dims <- paste0(query_reduction, 1:2)
  
  ref_df <- get_plot_df(ref, ref_reduction, ref_dims)
  query_df <- get_plot_df(query, query_reduction, query_dims)
  
  plot_expr <- function(df, dims) {
    xvar <- sym(dims[1])
    yvar <- sym(dims[2])
    
    base <- ggplot(df, aes(x = !!xvar, y = !!yvar)) +
      theme_minimal(base_size = 9) +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.width = unit(0.5, "lines"),
        legend.spacing.x = unit(0.2, "lines"),
        legend.margin = margin(t = 1, b = 1)
      )
    
    if (grey_zero) {
      base <- base +
        geom_point(data = df[df$expr == 0, ], color = "lightgrey", size = 0.2) +
        geom_point(data = df[df$expr > 0, ],
                   aes(x = !!xvar, y = !!yvar, color = expr), size = 0.6)
    } else {
      base <- base + geom_point(aes(color = expr), size = 0.6)
    }
    
    base <- base + scale_color_viridis_c(option = viridis_option, name = "Expression")
    
    if (!is.null(xlim)) base <- base + xlim(xlim)
    if (!is.null(ylim)) base <- base + ylim(ylim)
    
    return(base)
  }
  
  p1 <- plot_expr(ref_df, ref_dims)
  p2 <- plot_expr(query_df, query_dims)
  
  final_plot <- layout_op(p1, p2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = goi,
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  
  dims <- if (layout == "vertical") c(7, 14) else c(14, 7)
  
  ggsave(
    filename = paste0(filename_prefix, cleanup(goi), ".png"),
    plot = final_plot,
    width = dims[1],
    height = dims[2]
  )
}

################################################################################
################################################################################
################################################################################



# check: is a seurat object layered or not?
check_layered_object <- function(seurat_obj) {
  assays <- names(seurat_obj@assays)  # Get assay names as character vector
  result <- lapply(assays, function(a) {
    assay_obj <- seurat_obj[[a]]
    if ("Layers" %in% slotNames(assay_obj)) {
      layers <- Layers(assay_obj)
      if (length(layers) > 1) {
        return(paste(a, "has layers:", paste(names(layers), collapse = ", ")))
      }
    }
    return(NULL)
  })
  
  result <- Filter(Negate(is.null), result)
  
  if (length(result) > 0) {
    cat("✅ Layered structure detected:\n")
    cat(paste(result, collapse = "\n"))
    return(TRUE)
  } else {
    cat("❌ No layered structure detected.\n")
    return(FALSE)
  }
}



################################################################################
################################################################################
################################################################################


#integration preflight - marginally useful
int_preflight <- function(object) {
  DefaultAssay(object)
  DefaultAssay(object)
  validObject(object)
  class(object)
  class(object[["RNA"]])
  class(object[["SCT"]])
  check_layered_object(object)
}


################################################################################
################################################################################
################################################################################



#robust data writer for making 3d umap movies
#this checks that you actually have a 3d reduction somewhere. 
# it defaults to looking for a umap with 3+ dimensions
# it allows you to pull some or all metadata along with coordinates
# it tells you what metadata was pulled
# if you fail to provide a seurat object with 3d reductions, it will offer to make the reductions for you with seurats
# RunUMAP function
# gives graceful errors


generate_3d_umap_data <- function(
    seurat_obj,
    output_file = "umap3d_plotdata.csv",
    metadata_cols = "seurat_clusters"
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a valid Seurat object.")
  }
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The 'Seurat' package is required but not installed.")
  }
  
  if (!"SCT" %in% Assays(seurat_obj)) {
    warning("'SCT' assay not found; using default assay.")
  } else {
    DefaultAssay(seurat_obj) <- "SCT"
  }
  
  # Check or prompt for 3D UMAP
  if (!"umap3d" %in% names(seurat_obj@reductions)) {
    if (interactive()) {
      answer <- readline(prompt = "3D UMAP embedding not found. Run it now? [y/n]: ")
      if (tolower(answer) == "y") {
        tryCatch({
          seurat_obj <- RunUMAP(
            seurat_obj,
            dims = 1:30,
            n.components = 3L,
            reduction.name = "umap3d"
          )
        }, error = function(e) {
          stop("Failed to run 3D UMAP: ", e$message)
        })
      } else {
        stop("3D UMAP embedding is required and was not found.")
      }
    } else {
      stop("3D UMAP embedding not found and cannot prompt to run in non-interactive mode.")
    }
  }
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  }
  
  # Validate metadata column selection
  if (length(metadata_cols) == 1 && metadata_cols == "all") {
    meta_cols <- colnames(seurat_obj@meta.data)
    message("Including all metadata columns: ", paste(meta_cols, collapse = ", "))
  } else {
    meta_cols <- metadata_cols
    missing <- setdiff(meta_cols, colnames(seurat_obj@meta.data))
    if (length(missing) > 0) {
      stop("The following metadata columns were not found in the Seurat object: ",
           paste(missing, collapse = ", "))
    }
  }
  
  # Ensure UMAP_1,2,3 exist
  umap_vars <- c("UMAP_1", "UMAP_2", "UMAP_3")
  if (!all(umap_vars %in% colnames(Embeddings(seurat_obj, "umap3d")))) {
    stop("3D UMAP reduction does not contain expected components UMAP_1, UMAP_2, UMAP_3")
  }
  
  # Extract data
  vars_to_fetch <- c(umap_vars, meta_cols)
  plot_data <- tryCatch({
    FetchData(seurat_obj, vars = vars_to_fetch)
  }, error = function(e) {
    stop("Failed to fetch data: ", e$message)
  })
  
  plot_data$label <- rownames(plot_data)
  
  tryCatch({
    write.csv(plot_data, file = output_file, row.names = FALSE)
  }, error = function(e) {
    stop("Failed to write CSV: ", e$message)
  })
  
  message("3D plot data written to ", output_file)
  return(invisible(plot_data))
}



################################################################################
################################################################################
################################################################################


#simple plotting function to highlight clusters and color by a samplesource
#fxn is specific to this analysis

make_cluster_highlight_plot <- function(obj, cluster_id, group_col = "snn_res.0.8", source_col = "samplesource") {
  obj <- SetIdent(obj, value = group_col)
  cells_in_cluster <- WhichCells(obj, idents = cluster_id)
  cell_sources <- obj[[source_col]][cells_in_cluster, 1]
  
  atlas_cells <- cells_in_cluster[cell_sources == "atlas"]
  culture_cells <- cells_in_cluster[cell_sources == "cultures"]
  
  DimPlot(
    obj,
    group.by = group_col,
    cells.highlight = list(
      Atlas = atlas_cells,
      Cultures = culture_cells
    ),
    cols.highlight = c("dodgerblue3", "firebrick2"),
    cols = "grey90"
  ) + ggtitle(paste("Cluster", cluster_id, "- Atlas (Blue), Cultures (Red)"))
}


# WarnPrompt(): minimal yes/no gate for big downloads etc


WarnPrompt <- function(
    message                   = "your warning message here",
    default_yes_noninteractive = TRUE,
    abort_message             = "Aborted by user."
) {
  default <- if (default_yes_noninteractive) TRUE else FALSE
  
  if (isFALSE(utils::askYesNo(message, default = default))) {
    stop(abort_message, call. = FALSE)
  }
  
  invisible(TRUE)   # return invisibly so it can be used as a guard
}



################################################################################
################################################################################
################################################################################



# Flip one or both axes of any stored dimensional reduction
flipper <- function(x, embedding = "umap", dims = NULL) {
  # Case 1: Seurat object — flip embeddings
  if (inherits(x, "Seurat")) {
    if (!embedding %in% names(x@reductions)) {
      stop("Reduction '", embedding, "' not found in the Seurat object.")
    }
    
    emb <- x[[embedding]]@cell.embeddings
    
    # Default: flip all dimensions
    if (is.null(dims)) {
      dims_idx <- seq_len(ncol(emb))
    } else {
      dims_idx <- if (is.character(dims)) {
        match(dims, colnames(emb), nomatch = 0L)
      } else {
        dims
      }
      if (any(dims_idx == 0L)) stop("Invalid dimension names in `dims`.")
    }
    
    # Flip
    emb[, dims_idx] <- -emb[, dims_idx]
    
    # Return flipped embedding matrix (does not modify object)
    return(emb)
  }
  
  # Case 2: ggplot — flip plot axes
  if (inherits(x, "gg")) {
    if (is.null(dims) || all(dims %in% c(1, 2, "x", "y"))) {
      flip_x <- any(dims %in% c(1, "x"))
      flip_y <- any(dims %in% c(2, "y"))
      
      if (is.null(dims)) flip_x <- flip_y <- TRUE
      
      if (flip_x && flip_y) return(x + scale_x_reverse() + scale_y_reverse())
      if (flip_x) return(x + scale_x_reverse())
      if (flip_y) return(x + scale_y_reverse())
      
      return(x)
    }
    stop("For ggplot input, `dims` must be 1, 2, 'x', 'y', or NULL.")
  }
  
  stop("`flipper()` only supports Seurat or ggplot objects.")
}



################################################################################
################################################################################
################################################################################



#helper to generate and write labeled cluster comparison
#modal transferred annotation per cluster
#helper to generate and write labeled modal cluster comparison
write_modal_summary <- function(seurat_obj, annotation_col, cluster_col, label) {
  #coerce to vectors in case column is a data frame
  annot_vec <- as.vector(seurat_obj[[annotation_col]])
  clust_vec <- as.vector(seurat_obj[[cluster_col]])
  
  #build contingency table
  tab <- table(annot_vec, clust_vec)
  tab_df <- as.data.frame.matrix(tab)
  
  #identify modal label and percentage for each cluster
  modal_id <- apply(tab_df, 2, function(col) rownames(tab_df)[which.max(col)])
  modal_pct <- apply(tab_df, 2, function(col) 100 * max(col) / sum(col))
  
  #build labeled summary
  modal_summary <- rbind(modal_id, sprintf("%.1f%%", modal_pct))
  rownames(modal_summary) <- c("modal_annotation", "modal_percent")
  
  #combine and write to file
  full_table <- rbind(modal_summary, as.matrix(tab_df))
  outfile <- here::here("tables", "annotation_comparisons",
                        paste0("cluster_by_", cluster_col, "_vs_", annotation_col, "_", label, ".tsv"))
  write.table(full_table, file = outfile, sep = "\t", quote = FALSE, col.names = NA)
  
  return(full_table)
}




################################################################################
################################################################################
################################################################################


##### plotting functions for making more query-ref plots
plot_allcultured_orig_ident <- function(output_dir = NULL, filename_prefix = "orig.ident") {
  outdir <- output_dir %||% here::here("plots", "QueryRefPlots", "DimPlots", "all_cultured", "orig.ident")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  atlas_df <- as.data.frame(Embeddings(atlas, "umap")[, 1:2])
  colnames(atlas_df) <- c("UMAP_1", "UMAP_2")
  
  cc_df <- as.data.frame(Embeddings(cc_mapped, "ref.umap")[, 1:2])
  cc_df$group <- cc_mapped$orig.ident
  colnames(cc_df)[1:2] <- c("UMAP_1", "UMAP_2")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.3) +
    geom_point(data = cc_df, aes(x = UMAP_1, y = UMAP_2, color = factor(group)), size = 0.5, alpha = 0.8) +
    scale_x_reverse() + scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(outdir, paste0(filename_prefix, ".png")), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(outdir, paste0(filename_prefix, "_with_legend.png")), p, width = 7, height = 7)
  
  legend <- cowplot::get_legend(p + theme(legend.position = "right"))
  ggsave(file.path(outdir, paste0(filename_prefix, "_legend.png")), cowplot::plot_grid(legend), width = 3, height = 3)
}


plot_allcultured_unint080 <- function(output_dir = NULL, filename_prefix = "unintegrated_cc_clusters0.8") {
  outdir <- output_dir %||% here::here("plots", "QueryRefPlots", "DimPlots", "all_cultured", filename_prefix)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  atlas_df <- as.data.frame(Embeddings(atlas, "umap")[, 1:2])
  colnames(atlas_df) <- c("UMAP_1", "UMAP_2")
  
  cc_df <- as.data.frame(Embeddings(cc_mapped, "ref.umap")[, 1:2])
  cc_df$group <- cc_mapped$unintegrated_cc_clusters0.8
  colnames(cc_df)[1:2] <- c("UMAP_1", "UMAP_2")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.3) +
    geom_point(data = cc_df, aes(x = UMAP_1, y = UMAP_2, color = factor(group)), size = 0.5, alpha = 0.8) +
    scale_x_reverse() + scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(outdir, paste0(filename_prefix, ".png")), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(outdir, paste0(filename_prefix, "_with_legend.png")), p, width = 7, height = 7)
  
  legend <- cowplot::get_legend(p + theme(legend.position = "right"))
  ggsave(file.path(outdir, paste0(filename_prefix, "_legend.png")), cowplot::plot_grid(legend), width = 3, height = 3)
}


plot_allcultured_unint1 <- function(output_dir = NULL, filename_prefix = "unintegrated_cc_clusters1") {
  outdir <- output_dir %||% here::here("plots", "QueryRefPlots", "DimPlots", "all_cultured", filename_prefix)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  atlas_df <- as.data.frame(Embeddings(atlas, "umap")[, 1:2])
  colnames(atlas_df) <- c("UMAP_1", "UMAP_2")
  
  cc_df <- as.data.frame(Embeddings(cc_mapped, "ref.umap")[, 1:2])
  cc_df$group <- cc_mapped$unintegrated_cc_clusters1
  colnames(cc_df)[1:2] <- c("UMAP_1", "UMAP_2")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.3) +
    geom_point(data = cc_df, aes(x = UMAP_1, y = UMAP_2, color = factor(group)), size = 0.5, alpha = 0.8) +
    scale_x_reverse() + scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(outdir, paste0(filename_prefix, ".png")), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(outdir, paste0(filename_prefix, "_with_legend.png")), p, width = 7, height = 7)
  
  legend <- cowplot::get_legend(p + theme(legend.position = "right"))
  ggsave(file.path(outdir, paste0(filename_prefix, "_legend.png")), cowplot::plot_grid(legend), width = 3, height = 3)
}


# predicted.devtime
plot_predicted_devtime <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.devtime")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.devtime.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.devtime`,
    alpha = score - 0.1536 #subtract minimum value
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.devtime.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.devtime_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.devtime_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_dev_states_AB <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.dev_states_AB")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.dev_states_AB.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.dev_states_AB`,
    alpha = score - 0.3355 #subtraact minimum value
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.dev_states_AB.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.dev_states_AB_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.dev_states_AB_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_celltype_AB <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.celltype_AB")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.celltype_AB.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.celltype_AB`,
    alpha = score - 0.2055 #subtract min score
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.celltype_AB.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.celltype_AB_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.celltype_AB_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_annotations_AJ <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.annotations_AJ")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  keep <- query_obj$`predicted.annotations_AJ` != "nodata"
  query_obj <- subset(query_obj, cells = colnames(query_obj)[keep])
  
  score <- query_obj$`predicted.annotations_AJ.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.annotations_AJ`,
    alpha = score - 0.1683 #subtract minimum value
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.annotations_AJ.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.annotations_AJ_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.annotations_AJ_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_unint1 <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.unintegrated_atlas_clusters1")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.unintegrated_atlas_clusters1.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.unintegrated_atlas_clusters1`,
    alpha = score - 0.3537 #subtract min score
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters1.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters1_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters1_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_unint080 <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.unintegrated_atlas_clusters0.8")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.unintegrated_atlas_clusters0.8.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.unintegrated_atlas_clusters0.8`,
    alpha = score - 0.3597 #subtract min score
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters0.8.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters0.8_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters0.8_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}


plot_predicted_dev_states_AB <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.dev_states_AB")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.dev_states_AB.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.dev_states_AB`,
    alpha = score^2
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.dev_states_AB.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.dev_states_AB_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.dev_states_AB_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_celltype_AB <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.celltype_AB")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.celltype_AB.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.celltype_AB`,
    alpha = score^2
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.celltype_AB.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.celltype_AB_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.celltype_AB_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_annotations_AJ <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.annotations_AJ")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  keep <- query_obj$`predicted.annotations_AJ` != "nodata"
  query_obj <- subset(query_obj, cells = colnames(query_obj)[keep])
  
  score <- query_obj$`predicted.annotations_AJ.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.annotations_AJ`,
    alpha = score^2
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.annotations_AJ.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.annotations_AJ_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.annotations_AJ_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_unint1 <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.unintegrated_atlas_clusters1")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.unintegrated_atlas_clusters1.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.unintegrated_atlas_clusters1`,
    alpha = score^2
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters1.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters1_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters1_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}

plot_predicted_unint080 <- function(query_obj, object_name = "cc_mapped") {
  dir_path <- here::here("plots", "QueryRefPlots", "DimPlots", object_name, "predicted.unintegrated_atlas_clusters0.8")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  score <- query_obj$`predicted.unintegrated_atlas_clusters0.8.score`
  df <- data.frame(
    x = Embeddings(query_obj, "ref.umap")[,1],
    y = Embeddings(query_obj, "ref.umap")[,2],
    group = query_obj$`predicted.unintegrated_atlas_clusters0.8`,
    alpha = score^2
  )
  
  atlas_df <- data.frame(Embeddings(atlas, "umap"))
  colnames(atlas_df) <- c("x", "y")
  
  p <- ggplot() +
    geom_point(data = atlas_df, aes(x = x, y = y), color = "grey80", size = 0.2) +
    geom_point(data = df, aes(x = x, y = y, color = group, alpha = alpha, size = 0.2)) +
    scale_alpha(range = c(0.1, 1), guide = "none") +
    scale_size(range = c(0.005, 3), guide = "none") +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_minimal()
  
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters0.8.png"), p + theme(legend.position = "none"), width = 7, height = 7)
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters0.8_with_legend.png"), p, width = 7, height = 7)
  legend <- cowplot::get_legend(p + guides(alpha = "none", size = "none") + theme(legend.position = "right"))
  ggsave(file.path(dir_path, "predicted.unintegrated_atlas_clusters0.8_legend.png"), cowplot::plot_grid(legend), width = 3, height = 3)
}




################################################################################
################################################################################
################################################################################




# Load necessary libraries (ensure these are loaded in your session)
# library(dplyr)
# library(ggplot2)
# library(here)
# library(svglite) # Recommended for high-quality SVG output

# --- Function Definition ---

#' Generate and Save Tiny Comparison Line Plot as SVG
#'
#' Creates a small line plot comparing cell counts over time ('devtime' for atlas,
#' rounded 'pseudotime' for cultures) for a specific cluster and saves it as SVG.
#'
#' @param seurat_obj A Seurat object containing the metadata.
#' @param cluster_id The numerical ID of the cluster to plot.
#' @param resolution_col The name of the metadata column containing the cluster IDs
#'                       (e.g., "snn_res.0.8").
#' @param save_dir The directory path where the SVG file should be saved.
#'                 The function will create the directory if it doesn't exist.
#'                 Use here::here() to construct this path robustly.
#' @param plot_width The width of the output SVG file in inches. Default is 1.
#' @param plot_height The height of the output SVG file in inches. Default is 0.18.
#'
#' @return Invisibly returns the full path to the saved SVG file, or NULL if
#'         plotting failed (e.g., no data for the cluster).
#' @examples
#' \dontrun{
#' # Ensure necessary libraries and objects (e.g., lvint_full) are loaded
#' library(here)
#' output_directory <- here("plots", "micro-linePlots")
#'
#' # Example for a single cluster
#' save_tiny_plot_svg(
#'   seurat_obj = lvint_full,
#'   cluster_id = 34,
#'   resolution_col = "snn_res.0.8",
#'   save_dir = output_directory
#' )
#'
#' # Example iterating over multiple clusters
#' cluster_range <- 0:3
#' plot_paths <- lapply(cluster_range, function(cl_id) {
#'   save_tiny_plot_svg(
#'     seurat_obj = lvint_full,
#'     cluster_id = cl_id,
#'     resolution_col = "snn_res.0.8",
#'     save_dir = output_directory
#'   )
#' })
#' }
save_tiny_plot_svg <- function(seurat_obj,
                               cluster_id,
                               resolution_col = "snn_res.0.8",
                               save_dir,
                               plot_width = 1,
                               plot_height = 0.18) {
  
  # --- 1. Input Validation and Data Preparation ---
  
  # Check if resolution column exists
  if (!resolution_col %in% colnames(seurat_obj@meta.data)) {
    stop("Specified resolution column '", resolution_col, "' not found in Seurat metadata.")
  }
  
  # Filter metadata for the specific cluster
  # Use !!sym() to correctly handle the column name passed as a string
  cluster_metadata <- seurat_obj@meta.data %>%
    filter(!!sym(resolution_col) == cluster_id)
  
  # Check if any cells were found for this cluster
  if (nrow(cluster_metadata) == 0) {
    warning("No cells found for cluster ", cluster_id, " in column '", resolution_col, "'. Skipping plot.")
    return(invisible(NULL))
  }
  
  # Process Atlas data (using devtime)
  cluster_atlas <- cluster_metadata %>%
    filter(samplesource == "atlas") %>%
    filter(!is.na(devtime)) %>%
    group_by(devtime) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    mutate(samplesource = "atlas", time_point = devtime) %>%
    select(time_point, cell_count, samplesource)
  
  # Process Cultures data (using rounded pseudotime)
  cluster_cultures_raw <- cluster_metadata %>%
    filter(samplesource == "cultures")
  
  # Check for and handle NA pseudotimes
  na_pseudotime_count <- sum(is.na(cluster_cultures_raw$pseudotime))
  if (na_pseudotime_count > 0) {
    warning("For cluster ", cluster_id, ": Filtered out ", na_pseudotime_count,
            " cells from 'cultures' source due to NA pseudotime.")
  }
  
  cluster_cultures <- cluster_cultures_raw %>%
    filter(!is.na(pseudotime)) %>%
    mutate(rounded_pseudotime = round(pseudotime)) %>%
    group_by(rounded_pseudotime) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    mutate(samplesource = "cultures", time_point = rounded_pseudotime) %>%
    select(time_point, cell_count, samplesource)
  
  # Combine the two datasets
  plot_data <- bind_rows(cluster_atlas, cluster_cultures)
  
  # Check if there's any data to plot after processing
  if (nrow(plot_data) == 0) {
    warning("No valid time points with cell counts found for cluster ", cluster_id,
            " after processing atlas and cultures data. Skipping plot.")
    return(invisible(NULL))
  }
  
  # Find overall max Y value for setting y-axis limit
  y_max <- max(plot_data$cell_count, 1, na.rm = TRUE) # Ensure y_max is at least 1
  
  # Define colors (Atlas=Green, Cultures=Red)
  source_colors_ggplot <- c(
    "atlas" = "#21908CFF",
    "cultures" = "#BB3754FF"
  )
  
  # --- 2. Create Plot ---
  
  p <- ggplot(plot_data, aes(x = time_point, y = cell_count, color = samplesource)) +
    geom_line(
      linewidth = 0.3 # EDIT line thickness for plot lines here (e.g., 0.3)
    ) +
    scale_color_manual(values = source_colors_ggplot) +
    scale_x_continuous(
      breaks = c(12, 24),        # Set specific major breaks
      minor_breaks = c(6, 18),   # Set specific minor breaks
      labels = c("", ""),         # Provide vector of empty strings matching length of breaks
      limits = c(0, 26)          # Set x-axis limits to 0-26
    ) +
    scale_y_continuous(
      breaks = y_max,            # Set tick only at the overall max value
      labels = y_max,            # Label the tick with the max value
      expand = expansion(mult = c(0, 0.1)), # Expand y-axis slightly at the top
      limits = c(0, NA) # Ensure y-axis starts at 0
    ) +
    coord_cartesian(clip = "off") + # Turn off clipping
    theme_minimal(base_size = 6) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.2),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.2, color = "black"),
      axis.ticks.length = unit(0.5, "mm"),
      axis.text.y = element_text(
        size = 5, angle = 0, hjust = 0.5, vjust = 1.0, color = "black"
      ),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.position = "none",
      plot.background = element_rect(fill = 'transparent', color = NA),
      panel.background = element_rect(fill = 'transparent', color = NA)
    )
  
  # --- 3. Save Plot ---
  
  # Create the save directory if it doesn't exist
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Construct filename programmatically
  filename <- paste0("c_", cluster_id, "_cellnumber_x_devtime_x_source_line.svg")
  full_save_path <- file.path(save_dir, filename) # Use file.path for robustness
  
  # Save the plot as SVG
  # Using svglite is recommended for better quality if available
  ggsave(
    filename = full_save_path,
    plot = p,
    device = svglite::svglite, # Use svglite device (ensure package is installed)
    # device = "svg", # Alternative basic SVG device
    width = plot_width,
    height = plot_height,
    units = "in",
    # dpi = 600, # DPI is less relevant for SVG but can be kept
    bg = "transparent"
  )
  
  # Return the save path invisibly
  invisible(full_save_path)
}

# --- Example Usage ---

# # Ensure necessary libraries and objects (e.g., lvint_full) are loaded
# library(here)
#
# # 1. Define the output directory using here()
# output_directory <- here::here("plots", "micro-linePlots")
#
# # 2. Define the range of clusters to iterate over
# cluster_range <- 0:3 # Example: clusters 0, 1, 2, 3
#
# # 3. Specify the resolution column name
# res_col_name <- "snn_res.0.8"
#
# # 4. (Optional but recommended) Create the main directory once before looping
# dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
#
# # 5. Loop through the clusters and generate plots
# # Using lapply to apply the function to each cluster ID in the range
# # It will return a list of the paths where plots were saved (or NULL if skipped)
# saved_plot_paths <- lapply(cluster_range, function(current_cluster_id) {
#   message("Generating plot for cluster: ", current_cluster_id) # Progress message
#   save_tiny_plot_svg(
#     seurat_obj = lvint_full, # Replace with your actual Seurat object name
#     cluster_id = current_cluster_id,
#     resolution_col = res_col_name,
#     save_dir = output_directory
#   )
# })
#
# # Print the list of saved file paths (or NULLs for skipped plots)
# print(saved_plot_paths)

################################################################################
################################################################################
################################################################################



summarize_cluster_celltypes <- function(seurat_object,
                                        cluster_column_name,
                                        cluster_id,
                                        celltype_column_name = "predicted.celltype",
                                        score_column_name = "predicted.celltype.score") {

  meta_data <- seurat_object@meta.data

  required_cols <- c(cluster_column_name, celltype_column_name, score_column_name)
  missing_cols <- required_cols[!required_cols %in% colnames(meta_data)]
  
  # 4. Filter for Target Cluster
  target_cluster_metadata <- dplyr::filter(meta_data, .data[[cluster_column_name]] == cluster_id)

  # 5. Summary Table Generation
  summary_table <- dplyr::group_by(target_cluster_metadata, .data[[celltype_column_name]]) %>%
    dplyr::summarise(
      number_of_cells = dplyr::n(),
      mean_score = mean(.data[[score_column_name]], na.rm = TRUE),
      median_score = median(.data[[score_column_name]], na.rm = TRUE),
      q1_score = quantile(.data[[score_column_name]], 0.25, na.rm = TRUE),
      q3_score = quantile(.data[[score_column_name]], 0.75, na.rm = TRUE),
      min_score = min(.data[[score_column_name]], na.rm = TRUE),
      max_score = max(.data[[score_column_name]], na.rm = TRUE),
      sd_score = sd(.data[[score_column_name]], na.rm = TRUE),
    ) %>%
    dplyr::arrange(dplyr::desc(number_of_cells))

  # 6. Return Value
  return(summary_table)
}
