import os
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
                name="AbstractExtractor",
                func=lambda path: extract_abstract_from_pdf(path),
                description="useful for when you need to get the abstract from a scientific article pdf once i know the path to the pdf, which is pdf/(filename).pdf)",
            ),
            Tool(
                name="FilePathGetter",
                func=lambda i: get_sorted_pdf_file_path(i),
                description="useful for when you need to get the path to a pdf file in the pdf/ directory once i know the index of the file in the list of files in the pdf/ directory",
            ),
            Tool(
                name="Read Abstract",
                func=lambda path: extract_abstract_from_pdf(path),
                description="useful for when you need to get the abstract from a scientific article pdf once i know the path to the pdf, which is pdf/(filename).pdf)",
            ),
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
            self.tools, self.llm1, agent=AgentType.ZERO_SHOT_REACT_DESCRIPTION, verbose=True
        )
        self.worker = initialize_agent(
            self.tools, self.llm, agent=AgentType.ZERO_SHOT_REACT_DESCRIPTION, verbose=True
        )
    
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

    # Hardcode your question here if running wtih docker-compose
    #question = "Evaluation of effects of microplastics on zooplankton"

    # interactive input, not suitable for docker-compose run
    question = input("Enter your question: ")

    #question = "Effect of microplastics on different types of zooplankton?"
    files_with_answers = []
    # For every file in the pdf/ directory, extract the abstract and check if it contains the answer to the question, i is the number of all files in pdf dir
    for i in range(len(os.listdir('pdf/'))):  # Update directory path to 'pdf/'
        status = abstracter.search_abstracts(question, i)
        if status == 'True':
            files_with_answers.append(get_sorted_pdf_file_path(i))
        else:
            continue



    # Print the names of the PDFs that contain the possible answer in the abstract
    for file in files_with_answers:
        #Write into .txt file in output/ directory if it exists otherwise create it
        with open('output/output.txt', 'a') as f:
            f.write(file + '\n')
