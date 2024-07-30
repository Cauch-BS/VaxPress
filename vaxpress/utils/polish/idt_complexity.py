# Part of this code is identical to the SBOL-Utilities package
# This code was accesed at https://github.com/SynBioDex/SBOL-utilities?tab=License-1-ov-file
# This code was accessed on 2024-07-15 at commit 9331964
# The license for this code is the MIT License and is written below:
# MIT License

# Copyright (c) 2021 Raytheon BBN Technologies

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import asyncio
import aiohttp
import threading
from time import time
from requests import post
from requests.auth import HTTPBasicAuth
from typing import Dict, List

thread_lock = threading.Lock()

TOKEN_GEN_TIME = 0
ITERATION = 0
TOKEN = ''

class RateLimiter:
    def __init__(self, max_requests: int, time_interval: int) -> None:
        self.max_requests = max_requests
        self.time_interval = time_interval
        self.semaphore = asyncio.Semaphore(max_requests)
        self.calls = []

    async def acquire(self):
        await self.semaphore.acquire()
        current_time = time.time()
        if self.calls and current_time - self.calls[0] > self.time_interval:
            self.calls = [call for call in self.calls if current_time - call < self.time_interval]
        self.calls.append(current_time)
    
    def release(self):
        self.semaphore.release()

class IDTComplexity():
    '''Scoring function for complexity score from IDT.'''
    name = 'idt_complexity'
    description = 'Complexity Score from IDT'
    penalty_metric_flags = {'idt_complexity': 'i_comp'}
    
    _TOKEN_URL = "https://sg.idtdna.com/Identityserver/connect/token"
    _SCORE_URL = "https://sg.idtdna.com/restapi/v1/Complexities/ScreenGblockSequences"
    _TIMEOUT = 300 # Number of seconds to wait for score query requests to complete
    _MAX_LEN = 3000 # Maximum length of sequences to be considered by the IDT API
    _MIN_LEN = 125 # Minimum length of sequences to be considered by the IDT API
    _MAX_SEQ = 99 # Maximum number of sequences that can be accessed at once
    _MAX_REQ = 500 # Maximum number of requests that can be sent to the API server. 
    _N_WORKERS = 100 # Number of workers to send requests to the API server
    _MAX_RETRIES = 5 # Maximum number of retries for a request
    # _TOKEN_LIFE = 3200 # Number of seconds before the token expires

    def __init__(self, 
                #  weight: int = 1, interval : int = 1, _length_cds = None
                 ) -> None:
        # self.weight = weight 
        # self.interval = interval
        self.username = os.environ.get('IDT_USER', None)
        self.password = os.environ.get('IDT_PASSWD', None)
        self.client_id = os.environ.get('IDT_CLIENT', None)
        self.client_secret = os.environ.get('IDT_API', None)
        # self.priority = 37
        # global TOKEN_GEN_TIME, TOKEN, ITERATION
        # with thread_lock:
        #     self.iteration_count = ITERATION
        #     if int(time()) - TOKEN_GEN_TIME > self._TOKEN_LIFE:
        # self.token = self._get_access_token() 
        #     else: 
        #         self.token = TOKEN
    
    def _get_access_token(self) -> str:
        """Get access token for IDT API (see: https://www.idtdna.com/pages/tools/apidoc)
        :return: access token string
        """
        data = {
            'grant_type': 'password', 
            'username': self.username, 
            'password': self.password, 
            'scope': 'test'
            }
        auth = HTTPBasicAuth(self.client_id, self.client_secret)
        result = post(
            self._TOKEN_URL, data, auth=auth, timeout=self._TIMEOUT
        )

        if 'access_token' in result.json():
            global TOKEN_GEN_TIME, TOKEN
            with thread_lock:
                TOKEN_GEN_TIME = int(time())
                TOKEN = result.json()['access_token']
            return result.json()['access_token']
        raise ValueError("Error in response from IDT API. Please check credentials.")
    
    # @staticmethod
    # def parse_score(json_resp: List[List[Dict]]) -> List[int]:
    #     """Parse the JSON response from the IDT API to get the complexity score.
    #     :param json_resp: JSON response from IDT API"""
    #     scores = []
    #     for entries in json_resp:
    #         score = 0
    #         for entry in entries:
    #             if entry['IsViolated']:
    #                 score -= entry['Score']
    #             elif entry['Message'] == 'Authorization has been denied for this request.':
    #                 raise ValueError("Authorization has been denied for this request. Please check credentials.")
    #         scores.append(score)
    #     return scores
    
    async def _send_request(self, session, block: List[Dict]) -> List[Dict]:
        "Helper function to send an asynchronous API request to IDT."
        for attempt in range(self._MAX_RETRIES):
            try:
                async with session.post(self._SCORE_URL, json=block, timeout = self._TIMEOUT) as response:
                    response.raise_for_status()
                    return await response.json()
            except (aiohttp.ClientError, asyncio.TimeoutError) as error:
                if attempt < self._MAX_RETRIES - 1:
                    await asyncio.sleep(2 ** attempt)
                else:
                    raise error
    
    async def scoring(self, seqs: List[str]):
        """Calculate the complexity score for a given sequence from the IDT API.
        This system uses the gBlock API, which is intended for sequences between 125 and 3000bp in length.
        Sequences not within this range will throw an exception.
        A complexity score betwee 0 and 10 is synthesizable, while a score greater or equal to 10 is not synthesizable.
        The maximum number of sequences that can be accessed every minute is 500.
        If your population is greater than 500, an error will be thrown."""
        queries = []
        checked_length= False

        # Check if the sequence length is within the acceptable range
        for i, seq in enumerate(seqs):
            idx = seqs.find('AAAAAAAAAA')
            seq = seq[:idx]
            if not checked_length:
                if len(seq) not in range(self._MIN_LEN, self._MAX_LEN):
                    raise ValueError(f"Sequence length must be between {self._MIN_LEN} and {self._MAX_LEN}")
                checked_length = True
            #set up a list of queries as dictionaries
            queries.append(
                {"Name": f"seq_{i}", "Sequence": seq}
            )
        
        # Split the queries into blocks of sequences for more efficient requests. 
        block_query = [
            queries[i:i + self._MAX_SEQ] for i in range(0, len(queries), self._MAX_REQ)
        ]

        queue = asyncio.Queue()
        results = []
        rate_limiter = RateLimiter(self._MAX_REQ * 0.8, 60)

        
        async def worker(queue: asyncio.Queue, session: aiohttp.ClientSession, results) -> None:
            while True:
                block = await queue.get()
                await rate_limiter.acquire()
                try:
                    result = await self._send_request(session, block)
                    results.append(result)
                finally:
                    rate_limiter.release()
                    queue.task_done()

        async with aiohttp.ClientSession(
            headers = {"Authorization": f"Bearer {self.token}",
                    "Content-Type": "application/json; charset=utf-8"}, 
            connector = aiohttp.TCPConnector(limit = self._N_WORKERS))  as session:
            workers = [asyncio.create_task(worker(queue, session, results))
                    for _ in range(self._N_WORKERS)]
            for block in block_query:
                await queue.put(block)
            await queue.join()
        
        for worker in workers:
            worker.cancel()
        
        return results

    def score(self, seqs):
        """Asynchronously calcualtes the complexity score for a given sequence from the IDT API."""
        # if self.iteration_count % self.interval == 0:
        #     results = asyncio.run(self.scoring(seqs))
        #     global penalties, idt_complexity_score
        #     with thread_lock:
        #         penalties = self.parse_score(results[0])
        #         idt_complexity_score = [penalty * self.weight for penalty in penalties]
        
        # metrics = {'idt_complexity': penalties}
        # scores = {'idt_complexity': idt_complexity_score}

        # with thread_lock:
        #     global ITERATION
        #     ITERATION += 1

        # return scores, metrics

        results = asyncio.run(self.scoring(seqs))
        return results[0]

