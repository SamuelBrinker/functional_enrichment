# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import json

#from KBaseReport.KBaseReportClient import KBaseReport
from samuelbrinkerfunctional_enrichment_update.Utils.FunctionalEnrichmentUtil import FunctionalEnrichmentUtil
#END_HEADER


class samuelbrinkerfunctional_enrichment_update:
    '''
    Module Name:
    samuelbrinkerfunctional_enrichment_update

    Module Description:
    A KBase module: samuelbrinkerfunctional_enrichment_update
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com//samuelbrinkerfunctional_enrichment_update.git"
    GIT_COMMIT_HASH = "ebf49d4453db496b1a670c73e235627537ae03db"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        #END_CONSTRUCTOR
        pass

    def run_samuelbrinkerfunctional_enrichment_update(self, ctx, params):
        '''
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        '''
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_samuelbrinkerfunctional_enrichment_update
        print('--->\nRunning run_enrichment_update\nparams:')
        print(params, 'params')
        print(ctx, 'ctx')
        print(json.dumps(params, indent=1))

        for key, value in params.items():
                if isinstance(value, str):
                    params[key] = value.strip()


        fe1_runner = FunctionalEnrichmentUtil(self.config)
        returnVal=[]
        returnVal = fe1_runner.run_fe1(params)
        #END run_samuelbrinkerfunctional_enrichment_update

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_enrichment_update return value ' +
                             'returnVal is not type dict as required.')
        # return the results

        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
