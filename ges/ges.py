"""GES algorithm.

Imported from the Pcalg package.
Author: Diviyan Kalainathan

.. MIT License
..
.. Copyright (c) 2018 Diviyan Kalainathan
..
.. Permission is hereby granted, free of charge, to any person obtaining a copy
.. of this software and associated documentation files (the "Software"), to deal
.. in the Software without restriction, including without limitation the rights
.. to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
.. copies of the Software, and to permit persons to whom the Software is
.. furnished to do so, subject to the following conditions:
..
.. The above copyright notice and this permission notice shall be included in all
.. copies or substantial portions of the Software.
..
.. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
.. IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
.. FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
.. AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
.. LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
.. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
.. SOFTWARE.
"""
import os
import uuid
import warnings
import networkx as nx
from shutil import rmtree
from cdt.causality.graph.model import GraphModel
from pandas import read_csv
from cdt.utils.Settings import SETTINGS
from cdt.utils.R import RPackages, launch_R_script


def message_warning(msg, *a, **kwargs):
    """Ignore everything except the message."""
    return str(msg) + '\n'


warnings.formatwarning = message_warning


class GES_with_score(GraphModel):
    r"""GES algorithm **[R model]**.

    Args:
        lambda (float): threshold value for variable selection.

    Default Parameters:
       + FILE: '/tmp/cdt_GES/data.csv'
       + LAMBDA: str(0.001)
       + OUTPUT: '/tmp/cdt_GES/result.csv'
       + OUTPUT2: '/tmp/cdt_GES/result2.csv'

    """

    def __init__(self, lambda_=1):
        """Init the model and its available arguments."""
        if not RPackages.pcalg:
            raise ImportError("R Package pcalg is not available.")

        super(GES_with_score, self).__init__()
        self.arguments = {'{FOLDER}': '/tmp/cdt_GES/',
                          '{FILE_TRAIN}': 'train_data.csv',
                          '{FILE_VALID}': 'valid_data.csv',
                          '{LAMBDA}': str(1),
                          '{OUTPUT}': 'result.csv',
                          '{OUTPUT2}': 'result2.csv'}
        self.lambda_ = lambda_

    def create_graph_from_data(self, data, **kwargs):
        """Apply causal discovery on observational data using GES.

        Args:
            data (pandas.DataFrame): DataFrame containing the data

        Returns:
            networkx.DiGraph: Solution given by the GES algorithm.
        """
        # Building setup w/ arguments.
        self.arguments['{LAMBDA}'] = str(self.lambda_)
        results = self._run_ges(data)

        return nx.relabel_nodes(nx.DiGraph(results),
                                {idx: i for idx, i in enumerate(data.columns)})

    def get_score(self, train_data, valid_data, **kwargs):
        """Apply causal discovery on observational data using GES and return
        a training and a validation score (likelihood without penalization).

        Args:
            data (pandas.DataFrame): DataFrame containing the data

        Returns:
            Training score and validation score.
        """
        # Building setup w/ arguments.
        self.arguments['{LAMBDA}'] = str(self.lambda_)
        results = self._run_ges_with_score(train_data,  valid_data)
        dag = results[0]
        train_score = results[1][0][0]
        val_score = results[1][1][0]

        return nx.relabel_nodes(nx.DiGraph(dag),
                                {idx: i for idx, i in enumerate(train_data.columns)}), train_score, val_score


    def _run_ges_with_score(self, train_data, valid_data):
        """Setting up and running GES with all arguments."""
        id = str(uuid.uuid4())
        os.makedirs('/tmp/cdt_GES' + id + '/')
        self.arguments['{FOLDER}'] = '/tmp/cdt_GES' + id + '/'

        def retrieve_result():
            return read_csv('/tmp/cdt_GES' + id + '/result.csv',
                            delimiter=',').values, read_csv('/tmp/cdt_GES' + id + '/result2.csv', delimiter=',').values
        try:
            train_data.to_csv('/tmp/cdt_GES' + id + '/train_data.csv', header=False, index=False)
            valid_data.to_csv('/tmp/cdt_GES' + id + '/valid_data.csv', header=False, index=False)
            ges_result = launch_R_script("{}/ges_with_score.R".format(os.path.dirname(os.path.realpath(__file__))),
                                         self.arguments, output_function=retrieve_result)
        # Cleanup
        except Exception as e:
            rmtree('/tmp/cdt_GES' + id + '')
            raise e
        except KeyboardInterrupt:
            rmtree('/tmp/cdt_GES' + id + '/')
            raise KeyboardInterrupt
        rmtree('/tmp/cdt_GES' + id + '')
        return ges_result


    def _run_ges(self, data, fixedGaps=None, verbose=True):
        """Setting up and running GES with all arguments."""
        id = str(uuid.uuid4())
        os.makedirs('/tmp/cdt_GES' + id + '/')
        self.arguments['{FOLDER}'] = '/tmp/cdt_GES' + id + '/'

        def retrieve_result():
            return read_csv('/tmp/cdt_GES' + id + '/result.csv', delimiter=',').values

        try:
            data.to_csv('/tmp/cdt_GES' + id + '/data.csv', header=False, index=False)
            ges_result = launch_R_script("{}/ges_with_score.R".format(os.path.dirname(os.path.realpath(__file__))),
                                         self.arguments, output_function=retrieve_result, verbose=verbose)
        # Cleanup
        except Exception as e:
            rmtree('/tmp/cdt_GES' + id + '')
            raise e
        except KeyboardInterrupt:
            rmtree('/tmp/cdt_GES' + id + '/')
            raise KeyboardInterrupt
        rmtree('/tmp/cdt_GES' + id + '')
        return ges_result
