# Abstracter

## Project Setup

This project analyzes abstracts of scientific papers and determines their relevance to a user-defined question. It utilizes OpenAI's language model for this purpose.

There are two main ways to use this setup:

- Interactive: Running it interactively where you will be prompted for input.
- ~~Docker-compose: Running it using docker-compose where you will need to hardcode your question into `main.py`.~~ #DEPRECATED 

## Prerequisites

- Python 3.8+ 
- Docker (only needed if using docker-compose)
- An OpenAI account with an API key. If you don't have one, you can create it [here](https://platform.openai.com/account/api-keys)
- pubmed account
- Scopus API key
- IEEE API key

## Setup Instructions

1. Clone this repository:

`git clone git@github.com:grapergrape/Abstracter.git`


2. Open project using VSCode


3. Rename the `.envtemplate` file to `.env`:


4. Open the `.env` file and replace `sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx` with your OpenAI API key and add your pubmed, scopus and IEEE account/api key aswell.

## Running the Project

### Interactive Mode

Install the required dependencies:

`pip install -r requirements.txt`


To run the project interactively, ensure you are in the root directory of the project. Run `python main.py` and when prompted, you have thre supproted modes: 
- mode 1: is the same as the original version that reads the abstract and anwsers your question based on the abstract of the articles in question,
- mode 2: you hardcode the intial query that will use PubMeds and IEEE algorithms for querying relevant articles, and a mroe specific question that your LLM agent will use to go though abstracts of fetched articles, the output is stored in the corresponding JSON file,
- mode 3: you ask the LLM model a specific question about an article, he reads the scientific paper in chunks extracts relevant info from chunks, combines the chunk summaries into a coherent anwser based on the initial question> i got a detailed resposne from a 37 page scientific article in 1 minute and 20 seconds.

### Docker-compose Mode DEPRICATED due to interactive mode selection
# you can still use this if you hardcode everything 
To run the project using docker-compose:

1. Install Docker and docker-compose.

2. In `main.py`, replace the line where the `question` variable is defined with your predefined question:


`question = "Effect of microplastics on different types of zooplankton?"`

3. Build and run the Docker container:

`docker-compose up --build -d`

## Model Selection

This project supports two models:

    gpt-3.5-turbo-0301: A cheaper option but may yield lower quality results.
    gpt-4: Provides better performance but is more expensive.

You can switch between them in the code by using llm or llm1 when calling the model.

To set your preferred model, navigate to the Abstracter class in main.py and modify the model_name parameter in the ChatOpenAI and OpenAI instances under the __init__ method.
If you want to use custom model see what models are supproted by LangChain

