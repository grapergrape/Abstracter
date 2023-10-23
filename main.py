import os
import re
from dotenv import load_dotenv
from langchain.agents import Tool, initialize_agent, AgentType
from langchain.agents import initialize_agent, AgentType
from langchain.memory import ConversationBufferMemory
from helpers.toolset import extract_abstract_from_pdf, get_sorted_pdf_file_path, store_text, read_text, inverse_boolean_string
from langchain.chat_models import ChatOpenAI
# Load environment variables from .env
load_dotenv()

# Abstracter class definition
class Abstracter:
    def __init__(self):
        self.memory = ConversationBufferMemory(memory_key="chat_history", return_messages=True)
        self.llm1 = ChatOpenAI(model_name="gpt-4")
        self.llm = ChatOpenAI(model_name="gpt-3.5-turbo-0301")
        self.tools = [
            Tool(
                name="StoreText",
                func=lambda string: store_text(string),
                description="useful for when you need to store text of any length to helperdir/helper.txt",
            ),
            Tool(
                name="Read Text",
                func=lambda dummy: read_text(dummy),
                description="useful for when you need to read from helperdir/helper.txt",
            )

        ]
        self.agent = initialize_agent(
            self.tools, self.llm1, agent=AgentType.ZERO_SHOT_REACT_DESCRIPTION, verbose=True, handle_parsing_errors=True
        )
        self.worker = initialize_agent(
            self.tools, self.llm, agent=AgentType.ZERO_SHOT_REACT_DESCRIPTION, verbose=True, handle_parsing_errors=True
        )
    
    def get_relevance_score(self, question, i):
        path = get_sorted_pdf_file_path(i)
        # Check if the path is valid
        if not os.path.exists(path):
            print("Error: Invalid file path.")
            return
        text = extract_abstract_from_pdf(path)

        # prompt the worker to assign a relevance score to the text based on question
        relevance_prompt = f"If the text stored in helperdir/helper.txt is entirely relevant to the question '{question}', it would be 100%. If it is not relevant at all, it would be 0%. Make up and output a number to the relevance of the text to the question in percentages?"
        relevance_output = self.worker.run(relevance_prompt)
        relevance_score = int(re.search(r'\d+', relevance_output).group())
        print()
        # prompt the worker to write an explanation for the relevance score
        explanation_prompt = f"Explain, in a single paragraph, why you assigned the relevance score of {relevance_score} to the text: stored in helperdir/helper.txt based on the question: {question}."
        explanation = self.worker.run(explanation_prompt)

        return relevance_score, explanation
    
    def search_abstracts(self, question, i):
        path = get_sorted_pdf_file_path(i)
        # Check if the path is valid
        if not os.path.exists(path):
            print("Error: Invalid file path.")
            return
        text = extract_abstract_from_pdf(path)

        # Change between worker (cheapest gpt 3) and agent (gpt4)
        status_prompt = f"Analyze the text stored in helperdir/helper.txt. Does it mentiion anything that would indicate a relation to this: (True/False): {question}?"
        status = self.agent.run(status_prompt)


        # for explanation on decision making uncomment this
        """reasoning_prompt = f"Explain why you assigned the value of {status} for the relevance to the field of study: {question} on the text in helperdir/helper.txt."
        reasoning = self.agent.run(reasoning_prompt)"""

        # clear helper.txt
        store_text("")

        #print(reasoning)
        return status
    
if __name__ == "__main__":

    # Create instance of Abstracter class
    abstracter = Abstracter()   

    #Clear output.txt file
    with open('output/output.txt', 'w') as f:
        f.write('')

    # Clear report.md file
    with open('output/report.md', 'w') as f:
        f.write('')

    # Hardcode your question here if running wtih docker-compose
    question = "Evaluation of effects of microplastics on zooplanktons ability to feed it self"

    # interactive input, not suitable for docker-compose run
    #question = input("Enter your question: ")

    #question = "Effect of microplastics on different types of zooplankton?"
    files_with_answers = []
    article_scores = dict()  # Dictionary to store the relevance scores of each article

    # For every file in the pdf/ directory, extract the abstract and check if it contains the answer to the question, i is the number of all files in pdf dir
    for i in range(len(os.listdir('pdf/'))):  # Update directory path to 'pdf/'
        status = abstracter.search_abstracts(question, i)
        score, explanation = abstracter.get_relevance_score(question, i)
        article_scores[get_sorted_pdf_file_path(i)] = (score, explanation)
        if status == 'True':
            files_with_answers.append(get_sorted_pdf_file_path(i))
        else:
            continue


    # Sort the dictionary by score
    sorted_articles = dict(sorted(article_scores.items(), key=lambda item: item[1][0], reverse=True))




    # Print the names of the PDFs that contain the possible answer in the abstract
    for file in files_with_answers:
        #Write into .txt file in output/ directory if it exists otherwise create it
        with open('output/output.txt', 'a') as f:
            f.write(os.path.basename(file) + '\n')

    # Write the report
    with open('output/report.md', 'a') as f:
        first_below_80 = True
        for file, (score, explanation) in sorted_articles.items():
            if first_below_80 and score < 80:
                f.write('---\n')
                f.write('## Articles with Score Below 80\n')
                f.write('---\n')
                first_below_80 = False
            f.write(f'### {os.path.basename(file)}\n')
            f.write(f'Relevance Score: {score}\n')
            f.write(f'Explanation: {explanation}\n\n')