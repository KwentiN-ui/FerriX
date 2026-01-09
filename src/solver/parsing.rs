//! Helper Functions that make parsing INP-Input files easier.

/// This preprocess includes:
/// - removing leading and trailing whitespaces
/// - removing comments
/// - making all text uppercase
/// - merging lines that belong together (`,` at the end of line)
/// - removes empty lines
pub fn preprocess_inp(input_file: &str) -> String {
    input_file
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(|line| line.chars().map(|c| c.to_uppercase().to_string()).collect())
        // merge lines that end with `,`
        .map(|line: String| {
            if line.ends_with(',') {
                line + " "
            } else {
                line + "\n"
            }
        })
        // remove comments
        .filter(|line| !line.starts_with("**"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_preprocess_inp() {
        let inp = "**comment\n word  \n \t*keyword\n123.4\n4, 5, 6,\n7, 8, 9";
        assert_eq!(
            preprocess_inp(inp),
            "WORD\n*KEYWORD\n123.4\n4, 5, 6, 7, 8, 9\n"
        );
    }
}
