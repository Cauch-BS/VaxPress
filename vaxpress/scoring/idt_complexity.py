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

from . import ScoringFunction
from ..sequence import Sequence
from requests import post
from requests import HTTPBasicAuth
import asyncio
import aiohttp

import os

class IDTComplexity(ScoringFunction):
    name = 'idt_complexity'
    description = 'Complexity Score from IDT'

    arguments = [
        ('weight',
         dict(type=float, default=1.0, metavar='WEIGHT',
              help='scoring weight for remote repeats (default: 1.0)'))
    ]

    penalty_metric_flags = {'idt_complexity': 'i_comp'}
    
    _TOKEN_URL = "https://sg.idtdna.com/Identityserver/connect/token"
    _SCORE_URL = "https://sg.idtdna.com/restapi/v1/Complexities/ScreenGblockSequences"
    _TIMEOUT = 120 # Number of seconds to wait for score query requests to complete
    _MAX_LEN = 3000 # Maximum length of sequences to be considered by the IDT API
    _MIN_LEN = 125 # Minimum length of sequences to be considered by the IDT API
    _MAX_REQ = 99 # Maximum number of sequences that can be accessed at once

    def __init__(self, weight, min_length, _length_cds):
        self.weight = weight 
        self.min_length = min_length
        self.username = os.environ.get('IDT_USER', None)
        self.password = os.environ.get('IDT_PASSWD', None)
        self.client_id = os.environ.get('IDT_CLIENT', 'None')
        self.client_secret = os.environ.get('IDT_API', None)
        self.priority = 37
        self.token = self._get_access_token()
    
    
    def _get_idt_access_token(self) -> str:
        """Get access token for IDT API (see: https://www.idtdna.com/pages/tools/apidoc)
        :return: access token string
        """
        data = {
            'grant_type': 'password', 
            'username': self.username, 
            'password': self.password, 
            'scope': 'test'}
        auth = HTTPBasicAuth(self.client_id, self.client_secret)
        result = post(
            self._TOKEN_URL, data, auth=auth, timeout=self._TIMEOUT
        )

        if 'access_token' in result.json():
            return result.json()['access_token']
        else:
            raise Exception('''Access token for IDT API could not be generated. 
                             Check your credentials.
                             Make sure you exported your credentials as environmental variables:
                                IDT_USER = {your IDT username}
                                IDT_PASSWD = {your IDT password}
                                IDT_CLIENT = {your IDT client ID}
                                IDT_API = {your IDT API client password}''')
    
    @staticmethod
    def parse_score(json_resp):
        scores = []
        for entries in json_resp:
            score = 0
            for entry in entries:
                if entry['IsViolated']:
                    score -= entry['Score']
            scores.append(score)
        return scores
    
    async def _send_request(self, session, block):
        "Helper function to send an asynchronous API request to IDT."
        async with session.post(self._SCORE_URL, json =  block) as resp:
            resp_list = await resp.json()
            if len(resp_list) != len(block):
                raise Exception("Error in response from IDT API. Please try again.")
            return resp_list

    async def score(self, seqs):
        """Calculate the complexity score for a given sequence from the IDT API.
        This system uses the gBlock API, which is intended for sequences between 125 and 3000bp in length.
        Sequences not within this range will throw an exception.
        A complexity score betwee 0 and 10 is synthesizable, while a score greater or equal to 10 is not synthesizable.
        The maximum number of sequences that can be accessed every minute is 500.
        If your population is greater than 500, an error will be thrown."""
        penalties = []
        queries = []
        checked_length=False

        # Check if the sequence length is within the acceptable range
        for i, seq in enumerate(seqs):
            seq = seq[:- Sequence(seq).get_polyA()]
            if not checked_length:
                if len(seq) not in range(self._MIN_LEN, self._MAX_LEN):
                    raise Exception(f"Sequence length must be between {self._MIN_LEN} and {self._MAX_LEN}")
                checked_length = True
            #set up a list of queries as dictionaries
            queries.append(
                {"Name": f"seq {i}", "Sequece": seq}
            )
        
        # Split the queries into blocks of sequences for more efficient requests. 
        block_query = [
            queries[i:i + self._MAX_REQ] for i in range(0, len(queries), self._MAX_REQ)
        ]
        limit = 500 // self._MAX_REQ # Maximum number of sequences that can be accessed every minute is 500.

        async with aiohttp.ClientSession(
            headers = {"Authorization": f"Bearer {self.token}", 
                       "Content-Type": "application/json; charset=utf-8"},
                       timeout=self._TIMEOUT) as session:
            tasks = []
            for block in block_query:
                tasks.append(self._send_request(session, block))
                if len(tasks) > limit:
                    await asyncio.sleep(60)
            
            results = await asyncio.gather(*tasks)

        penalties = self.parse_score(results)
        idt_complexity_score = [penalty * self.weight for penalty in penalties]

        metrics = {'idt_complexity': penalties}
        scores = {'idt_complexity': idt_complexity_score}

        return scores, metrics

