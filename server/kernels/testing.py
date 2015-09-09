import unittest, json, sys, logging
from request.manager import TaskRequest
from manager import kernelMgr
from modules.configuration import MERRA_TEST_VARIABLES
from modules.utilities import wpsLog
wpsLog.addHandler( logging.StreamHandler(sys.stdout) )
wpsLog.setLevel(logging.DEBUG)

class KernelTests(unittest.TestCase):

    def setUp(self):
        self.test_point = [ -137.0, 35.0, 85000 ]
        self.test_time = '2010-01-16T12:00:00'
        self.operations = [ "time.departures(v0,t)", "time.climatology(v0,t,annualcycle)", "time.value(v0)" ]
        self.def_task_args =  { 'region': self.getRegion(), 'data': self.getData() }

    def tearDown(self):
        pass

    def getRegion(self, ipt=0 ):
        return '{"longitude": %.2f, "latitude": %.2f, "level": %.2f, "time":"%s" }' % (self.test_point[0]+5*ipt,self.test_point[1]-5*ipt,self.test_point[2],self.test_time)

    def getData(self, vars=[0]):
        var_list = ','.join( [ ( '"v%d:%s"' % ( ivar, MERRA_TEST_VARIABLES["vars"][ivar] ) ) for ivar in vars ] )
        data = '{"%s":[%s]}' % ( MERRA_TEST_VARIABLES["collection"], var_list )
        return data

    def getOp(self, op_index ):
        return [ self.operations[ op_index ] ]

    def getResultData(self, result_data ):
        if "error" in result_data:
            print>>sys.stderr, result_data["error"]
        return result_data.get('data',[])

    def getTaskArgs(self, op, ipt=0 ):
        task_args = { 'region': self.getRegion(ipt), 'data': self.getData() }
        task_args['operation'] = op
        return task_args

    def test01_cache(self):
        cache_level = 85000.0
        result = kernelMgr.run( TaskRequest( request={ 'region': { "level": cache_level }, 'data': self.getData() } ) )
        result_stats = result[0]['result']
        self.assertEqual( json.loads(result_stats['region']), { "lev": [cache_level] } )

    def test02_departures(self):
        test_result = [  -1.405364990234375, -1.258880615234375, 0.840728759765625, 2.891510009765625, -18.592864990234375,
                        -11.854583740234375, -3.212005615234375, -5.311614990234375, 5.332916259765625, -1.698333740234375,
                          8.750885009765625, 11.778228759765625, 12.852447509765625 ]
        task_args = self.getTaskArgs( op=self.getOp( 0 ) )
        result = kernelMgr.run( TaskRequest( request=task_args ) )
        result_data = self.getResultData(result[0])
        self.assertEqual( test_result, result_data[0:len(test_result)] )

    def test03_annual_cycle(self):
        test_result = [38.20164659288194, 40.60203721788194, 39.744038899739586, 37.738803439670136,
                       34.758260091145836, 32.372643364800346, 33.70814344618056, 35.980190700954864,
                       37.239708794487846, 38.93236626519097, 39.45347425672743, 35.83015611436632]
        task_args = self.getTaskArgs( self.getOp( 1 ), 1 )
        kernelMgr.persist()
        result = kernelMgr.run( TaskRequest( request=task_args ) )
        result_data = self.getResultData(result[0])
        self.assertEqual( test_result, result_data[0:len(test_result)] )

    def test04_value_retreval(self):
        test_result = 28.41796875
        task_args = self.getTaskArgs( self.getOp( 2 ), 2 )
        result = kernelMgr.run( TaskRequest( request=task_args ) )
        result_data = self.getResultData(result[0])
        self.assertEqual( test_result, result_data )

    def test05_multitask(self):
        test_results = [ [ -1.405364990234375, -1.258880615234375, 0.840728759765625 ], [ 48.07984754774306, 49.218166775173614, 49.36114501953125 ], 59.765625 ]
        task_args = self.getTaskArgs( op=self.operations )
        results = kernelMgr.run( TaskRequest( request=task_args ) )
        for ir, result in enumerate(results):
            result_data = self.getResultData(result)
            test_result = test_results[ir]
            if hasattr( test_result, '__iter__' ):  self.assertEqual( test_result, result_data[0:len(test_result)] )
            else:                                   self.assertEqual( test_result, result_data )