from wps import logger
from wps.processes import data_manager
from wps.processes import esgf_operation
from cdms2 import cdscan
import shlex
import uuid
from pywps import config
import os
import cdms2

class AggregateOperation(esgf_operation.ESGFOperation):
    def __init__(self):
        super(AggregateOperation, self).__init__()

    @property
    def title(self):
        return 'Aggregate'

    def __call__(self, data_manager, status):
        
        # First dodsrc thing just in case

        data_manager.metadata(self.input()[0])

        output_path = config.getConfigValue('server', 'outputPath', '/var/wps')
        new_file_name = '%s.nc' % (str(uuid.uuid4()),)
        new_file = os.path.join(output_path, new_file_name)
        
        fo = cdms2.open(new_file,"w")
        step = 20
        for v in self.input():
            status("LOOKING AT FILE: %s" % v.uri)
            f=cdms2.open(v.uri)
            V=f[v.var_name]
            sh = V.shape
            dtime = None
            kargs = {}
            final_units= None
            if self.domain is not None:
                status("WE HAVE A DOMAIN %s" % self.domain)
                # ok user passed a domain
                axes = V.getAxisList()
                for d in self.domain.dimensions:
                    isTime = False
                    for a in axes:
                        if d.name == "time" or (d.name == a.id and a.isTime()):
                            isTime = True
                            dtime = d
                            break
                    if not isTime:
                        if d.crs.name == "values":
                            if isinstance(d.start,basestring):
                                start = str(d.start)
                            else:
                                start = d.start
                            if isinstance(d.end,basestring):
                                end = str(d.end)
                            else:
                                end = d.end
                            kargs[d.name] = (start,end)
                        elif d.crs.name == "indices":
                            kargs[d.name] = slice(int(d.start),int(d.end))
                        else:
                            raise Exception("Unknown crs: %s" % d.crs.name)
            if dtime is not None:
                status("WILL READ PART OF THE FILE")
                t = V.getTime()
                t = t.clone()
                try:
                    if isinstance(d.start,basestring):
                        start = str(d.start)
                    else:
                        start = d.start
                    if isinstance(d.end,basestring):
                        end = str(d.end)
                    else:
                        end = d.end
                    start,end = t.mapInterval((start,end))
                except Exception,err:
                    # Time does not match anything in this file
                    status('MapInterval Failed %s %s with err %s' % (f,dtime,err))
                    continue
            else:
                status("WILL READ ALL THE FILE")
                start = 0
                end = len(V.getTime())
            for i in range(start,end,step):
                status("DUMPING %i to %i" % (i,i+step))
                status("KARGS: %s" % kargs)
                data = V(time=slice(i,i+step),**kargs)
                t = data.getTime()
                if final_units is None:
                    final_units = t.units

                else:
                    t.toRelativeTime(final_units)
                fo.write(data,id=v.var_name)
        fo.close()
        dap_url = self.create_dap_url(new_file_name)
        self.set_output(dap_url, v.var_name) 
