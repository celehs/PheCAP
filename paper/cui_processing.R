library(PheCAP)

# Perform majority voting on CUIs
# You may need to change "metamap_output/" to the actual path
# to the folder that contains the output from metamap
cui_list <- phecap_perform_majority_voting("metamap_output/")

# Extract CUI information from UMLS and write to a dictionary file
# You may need to change "dict_file.txt" to the location that
# you want to save the dictionary file for future use
# You may also need to change user, password, host, dbname
# based on your database connection
phecap_generate_dictionary_file(
  cui_list, dict_file = "dict_file.txt",
  user = "username", password = "password",
  host = "localhost", dbname = "umls")
