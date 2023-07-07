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

import requests
import threading
import csv
from queue import Queue
from feapder.network.user_agent import get
from bs4 import BeautifulSoup
from loguru import logger
from tenacity import retry, stop_after_attempt

import requests
import threading
import csv
from queue import Queue
from feapder.network.user_agent import get
from bs4 import BeautifulSoup
from loguru import logger
from tenacity import retry, stop_after_attempt

class DownloadFile(object):
    def __init__(self):
        self.url = "https://www.brenda-enzymes.org/ligand.php?brenda_ligand_id={}"
        self.headers = {
            'Host': 'www.brenda-enzymes.org',
            'Referer': 'https://www.brenda-enzymes.org/ligand.php?brenda_ligand_id=68',
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/112.0.0.0 Safari/537.36'
        }
        self.url_queue = Queue()
        self.save_queue = Queue()

    @retry(stop = stop_after_attempt(5))
    def query(self, url):
        self.headers['User-Agent'] = get()
        headers = self.headers
        response = requests.get(url, headers = headers, timeout = 30).content
        soup = BeautifulSoup(response, 'html.parser')
        h1_tag = soup.find('h1')
        h1_text = h1_tag.text if h1_tag else ""
        return h1_text

    def save(self):
        with open('Z02_All_Ligands_Names/results.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["ID", "H1 Text"])
            while True:
                save_item = self.save_queue.get()
                if save_item == 'DONE':
                    break
                writer.writerow([save_item['id'], save_item['content']])
                logger.success("Successfully Saved ID: " + save_item['id'])
                self.save_queue.task_done()

    def get_url(self):
        for i in range(100, 200):
            url = self.url.format(i)
            item = {"id": str(i), 'url': url}
            self.url_queue.put(item)
        self.url_queue.put('DONE')

    def get_content(self):
        while True:
            try:
                item = self.url_queue.get()
                if item == 'DONE':
                    self.save_queue.put('DONE')
                    break
                url = item['url']
                id = item['id']
                content = self.query(url)
                save_item = {'id': id, 'content': content}
                self.save_queue.put(save_item)
                self.url_queue.task_done()
            except Exception as e:
                logger.error("Error in get_content: " + str(e))

    def run(self):
        t_list = []
        t_get_url = threading.Thread(target = self.get_url)
        t_list.append(t_get_url)

        for i in range(1):
            t_get_content = threading.Thread(target = self.get_content)
            t_list.append(t_get_content)

        t_save = threading.Thread(target = self.save)
        t_list.append(t_save)

        for t in t_list:
            t.start()

        for q in [self.url_queue, self.save_queue]:
            q.join()


if __name__ == '__main__':
    down = DownloadFile()
    down.run()























