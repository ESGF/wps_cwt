import datetime
import os
import tempfile
import uuid

import cdutil
from pywps import config
import esgf

from wps import logger
from wps.processes import data_manager
from wps.processes import esgf_operation

class CDUtilOperation(esgf_operation.ESGFOperation):
    def __init__(self):
        super(CDUtilOperation, self).__init__()

    @property
    def title(self):
        return 'CDUtil Averager'

    def __call__(self, data_manager, status):

        domain = self.parameter('domain',required=False)

        status('Read input data')

        input_var = data_manager.read(self.input()[0])

        axes_arg = ''.join((str(input_var.getAxisIndex(x)) for x in axes.values))

        status('Averagering "%s" over dimensions "%s"' %
               (input_var.id,
                ', '.join(axes.values)))
        
        start = datetime.datetime.now() 

        new_var = cdutil.averager(input_var, axis=axes_arg)

        #status('Finished averaging in %s' % (datetime.datetime.now()-start,))

        new_var_name = '%s_avg_%s' % (input_var.id, '_'.join(axes.values))

        new_file_name = '%s.nc' % (str(uuid.uuid4()),)

        output_path = config.getConfigValue('server', 'outputPath', '/var/wps')

        new_file = os.path.join(output_path, new_file_name)

        #status('Writing output to "%s"' % (new_file,))

        data_manager.write('file://' + new_file, new_var, new_var_name)

        dap_url = self.create_dap_url(new_file_name)

        self.set_output(dap_url, new_var_name)
