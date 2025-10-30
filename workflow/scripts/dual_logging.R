# Dual logging utility for Snakemake R scripts
# Tees both stdout and message (stderr) to the log file AND keeps them on the console.

setup_dual_logging <- function(log_path) {
    if (!is.null(log_path) && length(log_path) > 0) {
        log_file <- log_path
        log_con <- file(log_file, open = "wt")
        
        # Create a custom message function that writes to both console and log
        assign("dual_message", function(...) {
            msg <- sprintf(...)
            # Write to console (stdout)
            cat(msg, "\n", file = stdout())
            # Write to log file
            cat(msg, "\n", file = log_con)
            flush(log_con)
        }, envir = .GlobalEnv)
        
        # Override the message function temporarily for this step
        assign("original_message", message, envir = .GlobalEnv)
        assign("message", dual_message, envir = .GlobalEnv)
        
        # Set warning level
        options(warn = 1)
        
        # Return cleanup function
        return(function() {
            # Restore original message function
            if (exists("original_message", envir = .GlobalEnv)) {
                assign("message", get("original_message", envir = .GlobalEnv), envir = .GlobalEnv)
                rm("original_message", envir = .GlobalEnv)
            }
            # Clean up dual_message function
            if (exists("dual_message", envir = .GlobalEnv)) {
                rm("dual_message", envir = .GlobalEnv)
            }
            close(log_con)
        })
    } else {
        # If no log file, just use regular console output
        return(function() {
            
        })
    }
}
