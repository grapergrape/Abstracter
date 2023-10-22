import re
import fitz
import os

def extract_abstract_from_pdf(pdf_file_path):
    # Open the PDF file
    with fitz.open(pdf_file_path) as pdf_file:
        # Initialize an empty string to hold the text of the PDF
        context = ""

        # Loop over the pages and extract the text
        for page_num in range(pdf_file.page_count):
            page = pdf_file[page_num]
            context += page.get_text("text")
        
        # List of patterns to search for
        patterns = ['Abstract', 'ABSTRACT', 'A B S T R A C T', 'abstract', 'Abstract.', 'a b s t r a c t']

        # Initialize match and found_pattern
        match = None
        found_pattern = None

        # Search for the patterns in the text
        for pattern in patterns:
            match = re.search(pattern, context)
            if match:
                # If a match is found, store the pattern and break the loop
                found_pattern = pattern.split() # split pattern into words
                break

        # If a match is found, extract the next 450 words
        if match:
            # Split the context into words
            words = context.split()
            
            # find the index of the first word of the found pattern
            try:
                pattern_index = words.index(found_pattern[0])
            except ValueError:
                return "Pattern not found in the list of words"
            
            # Take the next 450 words after the found pattern
            abstract_words = words[pattern_index + len(found_pattern) : pattern_index + len(found_pattern) + 450]
            
            # Join the words back into a string
            abstract = ' '.join(abstract_words)

            # Clean abstract text by stripping any leading/trailing whitespaces
            abstract = abstract.strip()

            # Store the abstract text (Assuming the function store_text() is defined somewhere)
            store_text(abstract)

            return abstract

        # If no pattern is found, return an appropriate message
        return "No abstract found in the provided PDF."
    
def get_sorted_pdf_file_path(i):
    # List all files in the "pdf" directory
    files = os.listdir('pdf/')
    
    # Filter the list to only include files that end with '.pdf'
    pdf_files = [file for file in files if file.endswith('.pdf')]
    
    # Sort the list in alphabetical order
    pdf_files.sort()
    # make sure i is an integer since llms pass strings in tools
    i = int(i)
    # Make sure the provided index 'i' is not out of the list bounds
    if i < 0 or i >= len(pdf_files):
        return "Error: Index out of bounds. Please provide an index between 0 and {}.".format(len(pdf_files) - 1)
    
    # Get the absolute path of the i-th file
    abs_path = os.path.abspath(os.path.join('pdf/', pdf_files[i]))
    
    # Return the absolute path
    return abs_path

def store_text(string):
    # this function stores text of any length to output/helper.txt
    with open("helperdir/helper.txt", "w") as f:
        f.write(string)

def read_text(dummy):
    # this function reads output/helper.txt
    with open("helperdir/helper.txt", "r") as f:
        return f.read()
    
def inverse_boolean_string(string):
    # this function inverts a boolean string if worker decieds to change their answer
    if string == "True":
        return "False"
    elif string == "False":
        return "True"
    else:
        return "Error: String is not a boolean."