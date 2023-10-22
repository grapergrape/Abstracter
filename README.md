# Abstracter

## Project Setup

This project analyzes abstracts of scientific papers and determines their relevance to a user-defined question. It utilizes OpenAI's language model for this purpose.

There are two main ways to use this setup:

- Interactive: Running it interactively where you will be prompted for input.
- Docker-compose: Running it using docker-compose where you will need to hardcode your question into `main.py`.

## Prerequisites

- Python 3.8+ 
- Docker (only needed if using docker-compose)
- An OpenAI account with an API key. If you don't have one, you can create it [here]([https://platform.openai.com/account/api-keys]).

## Setup Instructions

1. Clone this repository:

`git clone git@github.com:grapergrape/Abstracter.git`


2. Open project using VSCode


3. Rename the `.envtemplate` file to `.env`:


4. Open the `.env` file and replace `sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx` with your OpenAI API key.

## Running the Project

### Interactive Mode

Install the required dependencies:

`pip install -r requirements.txt`


To run the project interactively, ensure you are in the root directory of the project. Run `python main.py` and when prompted, enter your question.

### Docker-compose Mode

To run the project using docker-compose:

1. Install Docker and docker-compose.

2. In `main.py`, replace the line where the `question` variable is defined with your predefined question:


`question = "Effect of microplastics on different types of zooplankton?"`

3. Build and run the Docker container:

`docker-compose up --build -d`

## Model Selection

This project supports two models:

    text-davinci-003: A cheaper option but may yield lower quality results.
    gpt-4: Provides better performance but is more expensive.

To set your preferred model, navigate to the Abstracter class in main.py and modify the model_name parameter in the ChatOpenAI and OpenAI instances under the __init__ method.

## Output

The analysis results are written to output/output.txt. Each line of the file contains the file path to a PDF that was found to be relevant to the question.

