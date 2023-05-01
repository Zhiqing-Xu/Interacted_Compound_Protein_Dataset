#!/usr/bin/env python
# coding: utf-8
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
import requests 
import threading
#--------------------------------------------------#
# To manage the queue for threading
from queue import Queue
#--------------------------------------------------#
# To obtain user agent for requests
from feapder.network.user_agent import get
#--------------------------------------------------#
# For logging purposes
from loguru import logger
#--------------------------------------------------#
# For retry mechanism in case of failures
from tenacity import retry, stop_after_attempt


###################################################################################################################
###################################################################################################################

class DownloadFile(object):
    # Defines the constructor method for the Download File class, which initializes the object's properties.
    def __init__(self):
        self.url = "https://www.brenda-enzymes.org/molfile.php?LigandID={}"
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
        response = requests.get(url, headers=headers, timeout = 30).content
        return response

    def save(self):
        while True:
            save_item = self.save_queue.get()
            filename = save_item['id'] + ".mol"
            with open('./Z02_All_Ligands_Molfiles/' + filename, 'wb+') as fp:
                fp.write(save_item['content'])
            logger.success("Successfully Saved File: " + filename)
            self.save_queue.task_done()

    def get_url(self):
        for i in range(269999, 300000):
            url = self.url.format(i)
            item = {"id": str(i), 'url': url}
            self.url_queue.put(item)

    def get_content(self):
        while True:
            item = self.url_queue.get()
            url = item['url']
            id = item['id']
            content = self.query(url)
            save_item = {'id': id, 'content': content}
            self.save_queue.put(save_item)
            self.url_queue.task_done()

    def run(self):
        t_list = []
        t_get_url = threading.Thread(target = self.get_url)
        t_list.append(t_get_url)

        for i in range(5):
            t_get_content = threading.Thread(target = self.get_content)
            t_list.append(t_get_content)

        for i in range(5):
            t_save = threading.Thread(target = self.save)
            t_list.append(t_save)

        for t in t_list:
            t.start()

        for q in [self.url_queue, self.save_queue]:
            q.join()


if __name__ == '__main__':
    down = DownloadFile()
    down.run()
