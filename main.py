import os
import re
import json
from dotenv import load_dotenv

from langchain.agents import Tool, initialize_agent, AgentType
from langchain.memory import ConversationBufferMemory
from langchain.chat_models import ChatOpenAI

from helpers.toolset import extract_abstract_from_pdf, get_sorted_pdf_file_path, store_text, read_text, inverse_boolean_string
from helpers.get_abstracts import Fetcher
from helpers.results import get_results
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
    
    def get_relevance_score(self, question):
        # read helperdir/helper.txt and store it into text
        text = read_text("")
        # prompt the worker to assign a relevance score to the text based on question
        relevance_prompt = f'''If the text: {text} shows high possibility to conatin anwser to the question '{question}', it would be 100%. If it is not relevant at all, it would be 0%. 
                            Assign a relevance score to the text based on the question without using any tools, in a single step.
                            Always output an actual number (not X) to the relevance of the text to the question in percentages, which you will defend in the next prompt?'''
        relevance_output = self.worker.run(relevance_prompt)
        relevance_score = int(re.search(r'\d+', relevance_output).group())
        # prompt the worker to write an explanation for the relevance score
        explanation_prompt = f"Explain, in a single, detailed paragraph with comparison to the text {text} , why you assigned the relevance score of {relevance_score} to the text: {question}."
        explanation = self.worker.run(explanation_prompt)

        # clear helper.txt
        store_text("")

        return relevance_score, explanation
    
    def search_abstracts(self, question, i):
        path = get_sorted_pdf_file_path(i)
        # Check if the path is valid
        if not os.path.exists(path):
            print("Error: Invalid file path.")
            return
        text = extract_abstract_from_pdf(path)

        # Change between worker (cheapest gpt 3) and agent (gpt4)
        status_prompt = f"Analyze the text stored in helperdir/helper.txt. Does it indicate in anyway that it would have a relation to this: (True/False): {question}?"
        status = self.agent.run(status_prompt)


        # for explanation on decision making uncomment this
        """reasoning_prompt = f"Explain why you assigned the value of {status} for the relevance to the field of study: {question} on the text in helperdir/helper.txt."
        reasoning = self.agent.run(reasoning_prompt)"""


        #print(reasoning)
        return status
    
    def fetch_abstracts_online(self, query, nr):
        fetcher = Fetcher(query, nr)
        papers = fetcher.fetch()
        return papers
    
    def process_fetched_abstracts(self, papers, question):
        for paper in papers:
            title = paper['title']
            abstract = paper['abstract']
            store_text(abstract)
            status = self.agent.run(f"Analyze the text stored in helperdir/helper.txt. Does it mention anything that would indicate a relation to this: (True/False): {question}?")
            if status == "False":
                store_text("")
                continue
            else: 
                print(f"Title: {title}\nAbstract: {abstract}\nStatus: {status}\n\n\n")
            # if status true write into output/literature.txt the names of the articles
            if status == "True":
                with open('output/literature.txt', 'a') as f:
                    f.write(f"{title}\n")
            store_text("")

    def store_to_json(self, query, nr_results, papers):
        data = {
            "query": query,
            "number_of_results": nr_results,
            "papers": papers
        }
        with open('helperdir/pubmed.json', 'w') as f:
            json.dump(data, f, indent=2)

    def check_json(self, query, nr_results):
        # Check if pubmed.json exists
        if os.path.exists('pubmed.json'):
            # If it exists, load the data
            with open('pubmed.json', 'r') as f:
                data = json.load(f)
            # Compare the stored query and number of results to the provided parameters
            if data["query"] == query and data["number_of_results"] == nr_results and len(data["papers"]) > 0:
                # If they match and there are papers, return the papers
                return data["papers"]
        # If the file does not exist, or the stored data does not match the provided parameters, 
        # or there are no papers, return None
        return None

if __name__ == "__main__":

    # Create instance of Abstracter class
    abstracter = Abstracter()   
    #Clear output.txt file
    with open('output/output.txt', 'w') as f:
        f.write('')

    # Clear report.md file
    with open('output/report.md', 'w') as f:
        f.write('')

    mode = input("Enter:\n1 for local abstract analysis from pdf files, \n2 to fetch articles online, \n3 to get an anwser to a question from a pdf article. ")
    if mode == "1":
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
            score, explanation = abstracter.get_relevance_score(question)
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
                f.write(f"{os.path.basename(file)}\n")

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

    elif mode == "2":
        query = "Optical microplastic identification"
        # how many articles do you want to check from pubmd (note it searches in chronological order)
        nr_results = 20

        # Check if the same query has been rerun and OpenAI api failed after so we dont overwork pubmed api
        papers = abstracter.check_json(query, nr_results)

        if papers is None:
            papers = abstracter.fetch_abstracts_online(query, nr_results)
            # store the query and all papers in helperdir/pubmd.json for future use if openai api fails
            abstracter.store_to_json(query, nr_results, papers)

        if papers:
            question = "Performance analysis of different methods of optical classification of microplastics"
            abstracter.process_fetched_abstracts(papers, question)
    
    elif mode == "3":
        question = "What is flowcytometry and how does it preform?"
        anwser = get_results("pdf/as-74-9-1012-1.pdf", question)
        print(anwser)
        # Clear anwser.txt file
        with open('output/anwser.txt', 'w') as f:
            f.write('')
        # Write the anwser
        with open('output/anwser.txt', 'a') as f:
            f.write(f'{anwser}')
