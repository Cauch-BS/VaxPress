import os
import pprint
from requests import post
from requests.auth import HTTPBasicAuth
from typing import Dict, List


class IDTComplexity:
    """Scoring function for complexity score from IDT."""

    name = "idt_complexity"
    description = "Complexity Score from IDT"

    _TOKEN_URL = "https://sg.idtdna.com/Identityserver/connect/token"
    _SCORE_URL = "https://sg.idtdna.com/restapi/v1/Complexities/ScreenGblockSequences"
    _TIMEOUT = 300
    _MAX_LEN = 3000
    _MIN_LEN = 125
    _MAX_SEQ = 99

    def __init__(self):
        self.username = os.environ.get("IDT_USER", None)
        self.password = os.environ.get("IDT_PASSWD", None)
        self.client_id = os.environ.get("IDT_CLIENT", None)
        self.client_secret = os.environ.get("IDT_API", None)
        self.token = self._get_access_token()

    def _get_access_token(self) -> str:
        data = {
            "grant_type": "password",
            "username": self.username,
            "password": self.password,
            "scope": "test",
        }
        auth = HTTPBasicAuth(self.client_id, self.client_secret)
        result = post(self._TOKEN_URL, data=data, auth=auth, timeout=self._TIMEOUT)

        if "access_token" in result.json():
            return result.json()["access_token"]

        raise ValueError("Error in response from IDT API. Please check credentials.")

    def _send_request(self, block: List[Dict]) -> List[Dict]:
        headers = {
            "Authorization": f"Bearer {self.token}",
            "Content-Type": "application/json; charset=utf-8",
        }
        try:
            response = post(
                self._SCORE_URL, json=block, headers=headers, timeout=self._TIMEOUT
            )
            response.raise_for_status()
            return response.json()
        except Exception as error:
            raise error

    def scoring(self, seqs: List[str]):
        queries = []
        checked_length = False

        for i, seq in enumerate(seqs):
            idx = seq.find("AAAAAAAAAA")
            seq = seq[:idx]
            if not checked_length:
                if len(seq) not in range(self._MIN_LEN, self._MAX_LEN):
                    raise ValueError(
                        f"Sequence length must be between {self._MIN_LEN} and {self._MAX_LEN}. \n Current length is {len(seq)}."
                    )
                checked_length = True
            queries.append({"Name": f"seq_{i}", "Sequence": seq})

        block_query = [
            queries[i : i + self._MAX_SEQ]
            for i in range(0, len(queries), self._MAX_SEQ)
        ]

        results = []
        for block in block_query:
            result = self._send_request(block)
            results.append(result)

        return results

    def score(self, names: List[str], seqs: List[str]) -> dict[str, dict[str, object]]:
        analyzed = dict()
        scores = self.scoring(seqs)[0][:]
        for name, score in zip(names, scores):
            total_score = 0
            violated = dict()
            for dic in score:
                if dic["IsViolated"]:
                    violated[dic["Name"]] = {
                        "Message": dic["DisplayText"],
                        "Score": dic["Score"],
                    }
                    total_score += dic["Score"]
            formatted = pprint.pformat(violated, indent=4, width=100)
            analyzed[name] = {"Total Score": total_score, "Violated": formatted}
        return analyzed
