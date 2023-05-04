# SCC-alignment
repository for code associated with the SCC-aware alignment tutorial paper (Plaisier et al. 202x)


# Config file generation

Add your working directory path to the JSON file, e.g. using sed (note the use of "\" prior to symbols are necessary)

```
sed -i 's/\~/\/scratch\/working_dir/g' SCC-analysis_config.json
```

To check that the config is formatted properly, we suggest using a JSON validator tool, such as https://jsonlint.com/ 
