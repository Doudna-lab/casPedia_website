# CasPEDIA Internal Management
## Overview

<details>
<summary>CasPEDIA Framework</summary>
<ul>

CasPEDIA main parts are:
  <details>
    <summary>POSTGRES database</summary>
    <ul>
      
  - safely stores all the information loaded from the source information provided by the curators
  </details>
  
  <details>
    <summary>Python scripts</summary>
    <ul>
      
  - perform data parsing and formatting
  - these interact with the POSTGRES database to load and pull data as needed
  </details>
      
  <details>
    <summary>FLASK API</summary>
    <ul>
  
  - renders the web pages based on the data stored by the python scripts

  </details>
</details>

## Build and Update Data

<details>
<summary>Ingest CasPEDIA Master Table</summary>
<ul>

  - CasPEDIA relies on google forms spreadsheets hosted on google drive to load its Wiki pages
  - To capture those spreadsheets into the site API, we use three python scripts
  - The first, pulls and parses the `master table`, which contains highlevel content about the entries
```
python py/db_loadNupdate.py
```
</details>

<details>
<summary>Ingest CasPEDIA Entries</summary>
<ul>

  - Next we capture the entries' spreadsheets into the site API
  - In this step, we pull and parse the `individual entries`, which harbors more detailed content about the entries
  - This is where the CasPEDIA references are built into the database
```
python py/wiki_loadNupdate.py
```
</details>

<details>
<summary>Pre-load Information</summary>
<ul>
  
  - Next, the information parsed from the first step is stored into static pickled containers
  - These will be dynamically loaded onto HTML pages every time they're accessed through the website
  - By doing this, any changes to the page formatting (HTML) stays separated from the actual information (pickled content)
```
python py/preload_wiki.py
```
