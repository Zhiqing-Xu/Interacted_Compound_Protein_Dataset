#!/usr/bin/env python
# coding: utf-8
# author: Zhiqing Xu

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# The following code ensures the code work properly in 
# MS VS, MS VS CODE and jupyter notebook on both Linux and Windows.
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
if __name__ == "__main__":
    print("="*80)
    if os.name == 'nt' or platform == 'win32':
        print("Running on Windows")
        if 'ptvsd' in sys.modules:
            print("Running in Visual Studio")
#--------------------------------------------------#
    if os.name != 'nt' and platform != 'win32':
        print("Not Running on Windows")
#--------------------------------------------------#
    if "__file__" in globals().keys():
        print('CurrentDir: ', os.getcwd())
        try:
            os.chdir(os.path.dirname(__file__))
        except:
            print("Problems with navigating to the file dir.")
        print('CurrentDir: ', os.getcwd())
    else:
        print("Running in python jupyter notebook.")
        try:
            if not 'workbookDir' in globals():
                workbookDir = os.getcwd()
                print('workbookDir: ' + workbookDir)
                os.chdir(workbookDir)
        except:
            print("Problems with navigating to the workbook dir.")
#--------------------------------------------------#

###################################################################################################################
###################################################################################################################
# Imports
#--------------------------------------------------#
import re
import time
import copy
import pickle
import argparse
import numpy as np
import pandas as pd
#--------------------------------------------------#
import requests
import xmltodict
#--------------------------------------------------#
from timeit import timeit
#--------------------------------------------------#
import urllib
import xml.etree.ElementTree as ET
from urllib.request import urlopen

from bs4 import BeautifulSoup


###################################################################################################################
###################################################################################################################


from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.common.by import By
from webdriver_manager.firefox import GeckoDriverManager
from selenium.webdriver.firefox.options import Options

# Setup Firefox options
firefox_options = Options()

# Set your download path here
download_path = "C:/Users/cesiu/Desktop/github_all/Git_Interacted_Compound_Protein_Dataset/Z02_All_Ligands_Molfiles"

firefox_options.set_preference("browser.download.folderList", 2)
firefox_options.set_preference("browser.download.manager.showWhenStarting", False)
firefox_options.set_preference("browser.download.dir", download_path)
firefox_options.set_preference("browser.helperApps.neverAsk.saveToDisk", "chemical/x-mdl-molfile")

# Setup WebDriver
webdriver_service = Service(GeckoDriverManager().install())
driver = webdriver.Firefox(service=webdriver_service, options=firefox_options)

# URL of the webpage
url = 'https://www.brenda-enzymes.org/ligand.php?brenda_ligand_id=69' # replace with your URL

driver.get(url)

# Assuming the download link can be accessed by the class name 'download'
# find the download button and click it
download_button = driver.find_element(By.CLASS_NAME, 'download')
download_button.click()

# You should see a file being downloaded in the download_path directory
# Note: Selenium does not have built-in support to wait until the download is complete. 
# You may need to add a wait time or check the download directory for the completed download.

driver.quit()

























