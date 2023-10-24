import os
from Bio import Entrez
from dotenv import load_dotenv
import time

class Fetcher:
    def __init__(self, query, nr_results):
        load_dotenv()
        self.email = os.getenv('EMAIL')
        Entrez.email = self.email
        self.query = query
        self.nr_results = nr_results

    def search(self):
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax=self.nr_results,
                                retmode='xml',
                                term=self.query)
        results = Entrez.read(handle)
        time.sleep(0.1)  # Timeout to rate limit API key to not get IP blocked from NCBI
        return results

    def fetch_details(self, id_list):
        ids = ','.join(id_list)
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        results = Entrez.read(handle)
        time.sleep(0.1)  # Timeout to rate limit API key to not get IP blocked from NCBI
        return results

    def fetch(self):
        results = self.search()
        id_list = results['IdList']
        papers = self.fetch_details(id_list)
        
        papers_dict = []
        for i, paper in enumerate(papers['PubmedArticle']):
            title = paper['MedlineCitation']['Article']['ArticleTitle']
            try:
                abstract = paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            except KeyError:
                abstract = "No abstract available"
            papers_dict.append({"title": title, "abstract": abstract})
        
        return papers_dict


if __name__ == '__main__':
    fetcher = Fetcher('FTIR classification of microplastics', 5)
    papers = fetcher.fetch()
    print(papers)