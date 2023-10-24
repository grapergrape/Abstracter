import os
from Bio import Entrez
from dotenv import load_dotenv

class Fetcher:
    def __init__(self, query):
        load_dotenv()
        self.email = os.getenv('EMAIL')
        Entrez.email = self.email
        self.query = query

    def search(self):
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax='5',
                                retmode='xml',
                                term=self.query)
        results = Entrez.read(handle)
        return results

    def fetch_details(self, id_list):
        ids = ','.join(id_list)
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        results = Entrez.read(handle)
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
    fetcher = Fetcher('FTIR classification of microplastics')
    papers = fetcher.fetch()
    print(papers)
        
    