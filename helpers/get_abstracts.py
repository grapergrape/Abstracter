import os
import requests
from dotenv import load_dotenv
import time
import json

from Bio import Entrez

class Fetcher_pubmed:
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
            authors = [author['LastName'] + ' ' + author['ForeName'] for author in paper['MedlineCitation']['Article']['AuthorList']]
            try:
                abstract = paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            except KeyError:
                abstract = "No abstract available"
            papers_dict.append({"title": title, "abstract": abstract, "authors": authors})
        
        return papers_dict



class Fetcher_Scopus:

    def __init__(self, filename):
        load_dotenv()
        self.api_key = os.getenv('SCOPUS_API_KEY')
        try:
            with open(filename) as f:
                data = json.load(f)
            self.papers = data['papers']
            self.query = data.get('query',None)
            self.nr_results = data.get('number_of_results',None)
        except TypeError:
            self.papers = []
            self.query = None
            self.nr_results = None

    def search_paper(self, title):
        headers = {
            'Accept': 'application/json',
            'X-ELS-APIKey': self.api_key
        }
        params = {
            'query': 'TITLE("' + title + '")',
            'count': 1,
            'view': 'STANDARD'
        }
        response = requests.get('https://api.elsevier.com/content/search/scopus', headers=headers, params=params)
        time.sleep(0.1)
        return response.json().get('search-results', {}).get('entry', [])[0]

    def get_journal_metrics(self, source_title):
        headers = {
            'Accept': 'application/json',
            'X-ELS-APIKey': self.api_key
        }
        params = {
            'title': source_title
        }
        response = requests.get('https://api.elsevier.com/content/serial/title', headers=headers, params=params)
        time.sleep(0.1)
        
        journal = response.json().get('serial-metadata-response', {}).get('entry', [{}])[0]
        sjr = journal.get('SJRList', {}).get('SJR', [{}])[0].get('$', '')
        snip = journal.get('SNIPList', {}).get('SNIP', [{}])[0].get('$', '')
        citescore = journal.get('citeScoreYearInfoList', {}).get('citeScoreCurrentMetric', '')
        
        return sjr, snip, citescore

    def enrich_paper(self):
        for paper in self.papers:
            fetched_paper = self.search_paper(paper['title'])
            if fetched_paper == {}:
                print('Paper not found in Scopus')
                continue
            paper['citations_count'] = fetched_paper.get('citedby-count', [])
            paper['cover_date'] = fetched_paper.get('prism:coverDate', [])
            paper['source_title'] = fetched_paper.get('prism:publicationName', [])
            paper['impact_factor'] = self.get_journal_metrics(paper['source_title'])
            
            # Adding Scopus Id if not present
            scopus_id = fetched_paper.get('dc:identifier')
            if scopus_id:
                paper['scopus_id'] = scopus_id.split(':')[1]
            else:
                paper['scopus_id'] = []
        
        with open('helperdir/pubmed.json', 'w') as json_file:
            json.dump({"query": self.query, "number_of_results": self.nr_results,"papers": self.papers}, json_file, indent=4)
        
        return self.papers
    def h_index(self, author):
        # Split the author's name into first and last names
        first_name, last_name = author.split(' ')[:2]

        # Fetch author's papers
        headers = {
            'Accept': 'application/json',
            'X-ELS-APIKey': self.api_key
        }
        params = {
            'query': f'AUTHFIRST({first_name}) and AUTHLAST({last_name})',
            'count': 25,  # Adjust as needed
            'view': 'STANDARD'
        }

        response = requests.get('https://api.elsevier.com/content/search/author', headers=headers, params=params)
        time.sleep(0.1)  # Remember to respect the API rate limits
            
        if response.status_code != 200:
            print(f"Request failed with status {response.status_code}: {response.text}")
            return None

        author_papers = response.json().get('search-results', {}).get('entry', [])
        
        # Sort citation counts in descending order
        citation_counts = sorted([int(paper.get('citedby-count', 0)) for paper in author_papers], reverse=True)

        # Calculate h_score
        h_index = 0
        for i, citations in enumerate(citation_counts):
            if citations >= i + 1:
                h_index = i + 1
            else:
                break

        return h_index
class Fetcher_ieee:
    def __init__(self, query, nr_results):
        load_dotenv()
        self.api_key = os.getenv('IEEE_API_KEY')
        print(self.api_key)
        self.query = query
        self.nr_results = nr_results
        self.base_url = "http://ieeexploreapi.ieee.org/api/v1/search/articles"

    def search(self):
        params = {
            "apikey": self.api_key,
            "format": "json",
            "max_records": self.nr_results,
            "querytext": self.query
        }

        response = requests.get(self.base_url, params=params)
        print(response)
        results = response.json()

        # Time delay to adhere to rate limits
        time.sleep(0.1) 
        return results

    def fetch(self):
        results = self.search()
        papers_dict = []

        for paper in results['articles']:
            title = paper['title']
            authors = [author['full_name'] for author in paper['authors']]
            try:
                abstract = paper['abstract']
            except KeyError:
                abstract = "No abstract available"
            papers_dict.append({"title": title, "abstract": abstract, "authors": authors})

        return papers_dict
    
    def write_to_file(self, path='helperdir/ieee.json'):
        papers = self.fetch()
        output = {
            "query": self.query,
            "number_of_results": self.nr_results,
            "papers": papers
        }
        with open(path, 'w') as f:
            json.dump(output, f, ensure_ascii=False)


if __name__ == '__main__':
    """fetcher = Fetcher_ieee('microplastic optic detection', 15)
    fetcher.write_to_file()
    scopus_fetcher = Fetcher_Scopus('helperdir/ieee.json')
    enriched_papers = scopus_fetcher.enrich_paper()"""
    """fetcher = Fetcher_pubmed('FTIR classification of microplastics', 5)
    papers = fetcher.fetch()
    print(papers)"""
    h_score = Fetcher_Scopus('helperdir/pubmed.json').h_index('Jessica Caldwell')
    print(h_score)

