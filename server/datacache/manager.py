from modules.utilities import  *
from data_collections import getCollectionManger
from domains import *
from decomposition.manager import decompositionManager
from datacache.status_cache import StatusPickleMgr
from datacache.persistence.background_thread import PersistenceThread
import numpy, sys, os, traceback, Queue

import cdms2, time

def load_variable_region( dataset, name, cdms2_cache_args=None ):
    rv = None
    try:
        t0 = time.time()
        wpsLog.debug( "\n\n LLLLLLLOAD DataSet<%s> %x:%x, status = '%s', var=%s, args=%s \n\n" % ( dataset.id, id(dataset), os.getpid(), dataset._status_, name, str(cdms2_cache_args) ) )
        dset = dataset( name, **cdms2_cache_args )
        rv = numpy.ma.fix_invalid( dset )
        t1 = time.time()
        wpsLog.debug( " $$$ Variable '%s' %s loaded from dataset '%s' --> TIME: %.2f " %  ( name, str(rv.shape), dataset.id, (t1-t0) ) )
    except Exception, err:
        wpsLog.error( " ERROR loading Variable '%s' from dataset %s --> args: %s\n --- %s ----\n%s " %  ( name, str(dataset), str(cdms2_cache_args), str(err), traceback.format_exc() ) )
    return rv

def subset_variable_region( variable, cdms2_cache_args=None ):
    wpsLog.debug( " $$$ Subsetting Variable '%s' (%s): args='%s' " %  ( variable.id, str(variable.shape), str(cdms2_cache_args) ) )
    t0 = time.time()
    rv = numpy.ma.fix_invalid( variable( **cdms2_cache_args ) )
    t1 = time.time()
    wpsLog.debug( " >> Subset-> %s: TIME: %.2f " %  ( str(rv.shape), (t1-t0) ) )
    return rv

class CachedVariable:

    CACHE_NONE = 0
    CACHE_OP = 1
    CACHE_REGION = 2

    def __init__(self, **args ):
        self.id = args.get('id',None)
        self.type = args.get('type',None)
        self.dataset = args.get('dataset',None)
        self.specs = args
        self.stat = args.get('vstat',{})
        self.domainManager = DomainManager()

    def updateStats(self,tvar):
        self.stat['fill_value'] = tvar.fill_value
        self.stat['attributes'] = filter_attributes( tvar.attributes, [ 'units', 'long_name', 'standard_name', 'comment'] )
        self.stat['grid'] = tvar.getGrid()
        self.stat['id'] = tvar.id
        self.stat['cid'] = self.id
        self.stat['dtype'] = tvar.dtype
        cdms_domain = tvar.getDomain()
        self.stat['axes'] = [ d[0] for d in cdms_domain ]

    def getCacheSize(self, cache_map ):
        return self.domainManager.getCacheSize( cache_map )

    def persist( self, **args ):
        kwargs = { 'cid': self.id }
        kwargs.update( args )
        self.domainManager.persist( **kwargs )

    def stats( self, **args ):
        kwargs = { 'cid': self.id }
        kwargs.update( args )
        vstat = dict(self.stat)
        vstat['domains'] = self.domainManager.stats( **kwargs )
        return vstat

    @classmethod
    def getCacheType( cls, use_cache, operation ):
        if not use_cache: return cls.CACHE_NONE
        return cls.CACHE_REGION if operation is None else cls.CACHE_OP

    def addDomain(self, region, **args ):
        data = args.get( 'data', None )
        domain = Domain( region, tvar=data, vstat=self.stat )
        request_queue = args.get( 'queue', None )
        if request_queue: domain.cacheRequestSubmitted( request_queue )
        self.domainManager.addDomain( domain )
        return domain

    def addCachedDomain( self, domain_spec, **args ):
        domain = Domain( domain_spec, **args  )
        self.domainManager.addDomain( domain )
        return domain

    def uncache( self, region ):
        return self.domainManager.uncache( region )

    def cacheType(self):
        return self.specs.get( 'cache_type', None )

    def getSpec( self, name ):
        return self.specs.get( name, None )

    def findDomain( self, search_region  ):
        return self.domainManager.findDomain( Domain( search_region ) )

class CacheManager:
    RESULT = 1
    modifiers = { '@': RESULT };

    def __init__( self, name, **args ):
        self._cache = {}
        self.stat = {}
        self.name = name
        self.statusCache = StatusPickleMgr( '.'.join( [ 'stats_cache', name ] ) )
        self.loadStats()

    def persist( self, **args ):
        for cached_cvar in self._cache.values():
            cached_cvar.persist(**args)
        self.persistStats( **args )

    def uncache( self, datasets, region ):
        for dset in datasets:
            var_cache_id =  ":".join( [dset['collection'],dset['name']] )
            cached_cvar = self._cache.get( var_cache_id, None )
            if cached_cvar: cached_cvar.uncache( region )
            else: wpsLog.error( " Attempt to uncache unrecognized variable: %s" % var_cache_id )

    def stats( self, **args ):
        cache_stats = []
        for cached_cvar in self._cache.values():
            cache_stats.append( cached_cvar.stats(**args) )
        return cache_stats

    def loadStats( self, **args ):
        stats = self.statusCache['stats']
     #   wpsLog.debug( "\n ***-------> Load Stats: %s\n" % stats )
        if stats:
            for var_stats in stats:
                cache_id = var_stats.get('cid',None)
                cached_cvar = self._cache.get( cache_id, None )
                if cached_cvar is None:
                    cached_cvar = CachedVariable(  id=cache_id, vstat=var_stats )
                    self._cache[ cache_id ] = cached_cvar
                domain_stats = var_stats['domains']
                for domain_stat in domain_stats:
                    cached_cvar.addCachedDomain( domain_stat['region'], dstat=domain_stat, vstat=var_stats )

    def persistStats( self, **args ):
        stats = self.stats( **args )
     #   wpsLog.debug( "\n ***-------> Persist Stats[%s]:\n %s\n" % ( args.get('loc',""), stats ) )
        self.statusCache['stats'] = stats

    @classmethod
    def getModifiers( cls, variable_name ):
        var_type = cls.modifiers.get( variable_name[0], None )
        return ( var_type, variable_name if ( var_type == None ) else variable_name[1:] )

    def getVariable( self, cid, new_region ):
        cached_cvar = self._cache.get( cid, None )
        if cached_cvar is None: return Domain.DISJOINT, None
        status, domain = cached_cvar.findDomain( new_region )
        wpsLog.debug( "Searching for cache_id '%s', cache keys = %s, domain=%s, Found var = %s, Found domain = %s, status = %d" %
                      ( cid, str(self._cache.keys()), str(new_region), (cached_cvar is not None), (domain is not None), status ) )
        return status, domain

    def addVariable( self, cache_id, data, specs ):
        cached_cvar = self._cache.get( cache_id, None )
        if cached_cvar is None:
            var_type, var_id = self.getModifiers( cache_id )
            cached_cvar = CachedVariable( type=var_type, id=cache_id, specs=specs )
            cached_cvar.updateStats( data )
            self._cache[ cache_id ] = cached_cvar
        region = specs.get( 'region', {} )
        cached_cvar.addDomain( region, data=data )


    def getResults(self):
        return self.filterVariables( { 'type': CachedVariable.RESULT } )

    def filterVariables( self, filter ):
        accepted_list = [ ]
        for var_spec in self._cache.values():
            accepted = True
            for filter_item in filter.items():
                if ( len(filter) > 1 ) and (filter[0] is not None) and (filter[1] is not None):
                    spec_value = var_spec.getSpec( filter[0] )
                    if (spec_value <> None) and (spec_value <> filter[1]):
                        accepted = False
                        break
            if accepted: accepted_list.append( var_spec )

class DataManager:

    def __init__( self, name, **args ):
        self.cacheManager = CacheManager( name, **args )
        self.collectionManager = getCollectionManger( **args )
        self.persist_queue = Queue.Queue()
        self.persistenceThread = None
        enable_background_persist = args.get( 'background_persist', True )
        if enable_background_persist:
            self.persistenceThread = PersistenceThread( args=(self.persist_queue,) )
            self.persistenceThread.start()

    def close(self):
        self.collectionManager.close()

    def getName(self):
        return self.cacheManager.name

    def persist( self, **args ):
        self.cacheManager.persist( **args )

    def persistStats( self, **args ):
        self.cacheManager.persistStats( **args )

    def uncache( self, data, region ):
        self.cacheManager.uncache( data, region )

    def stats( self, **args ):
        return self.cacheManager.stats( **args )

    def loadVariable( self, data, region, cache_type ):
        data_specs = {}
        domain =  None
        dataset = None
        cache_region = region
        name =  data.get( 'name', None );
        variable = data.get('variable',None)
        t0 = time.time()
        wpsLog.debug( " #@@ DataManager:LoadVariable %s (time = %.2f), region = [%s], cache_type = %d" %  ( str( data ), t0, str(region), cache_type ) )
        if variable is None:
            collection = data.get('collection',None)
            if collection is not None:
                var_cache_id =  ":".join( [collection,name] )
                status, domain = self.cacheManager.getVariable( var_cache_id, region )
                if status is not Domain.CONTAINED:
                    dataset = self.loadFileFromCollection( collection, name )
                    data_specs['dataset'] = record_attributes( dataset, [ 'name', 'id', 'uri' ])
                else:
                    variable = domain.getVariable()
                    data_specs['dataset']  = domain.spec.get( 'dataset', None )
                    data_specs['region'] = domain.getRegion()
                data_specs['cid']  = var_cache_id
            else:
                wpsLog.debug( " $$$ Empty Data Request: '%s' ",  str( data ) )
                return None, data_specs

            if (variable is None):
                if (cache_type == CachedVariable.CACHE_NONE):
                    variable = dataset[name]
                    data_specs['region'] = region
                else:
                    load_region = decompositionManager.getNodeRegion( region ) if (cache_type == CachedVariable.CACHE_REGION) else region
                    cache_region = Region( load_region, axes=[ CDAxis.LATITUDE, CDAxis.LONGITUDE, CDAxis.LEVEL ] )
                    if dataset == None: dataset = self.loadFileFromCollection( collection, name )
                    variable = load_variable_region( dataset, name, cache_region.toCDMS() )
                    data_specs['region'] = cache_region
                    self.persist_queue.put( ( variable.data, data_specs ) )

            else:
                wpsLog.debug( " *** Loading data from cache, cached region= %s" % str(data_specs['region']) )
                if (cache_type == CachedVariable.CACHE_OP) and (region is not None):
                    variable = subset_variable_region( variable, region.toCDMS() )
                    data_specs['region'] = region

            data_specs['variable'] = record_attributes( variable, [ 'long_name', 'name', 'units' ] )
        else:
            data_specs['variable'] = record_attributes( variable, [ 'long_name', 'name', 'id', 'units' ]  )

        if (variable is not None) and (cache_type <> CachedVariable.CACHE_NONE):
            data_specs['cache_type'] = cache_type
            if (domain is None) or not domain.equals( cache_region):
                self.cacheManager.addVariable( var_cache_id, variable, data_specs )
        t1 = time.time()
        wpsLog.debug( " #@@ DataManager:FinishedLoadVariable %s (time = %.2f, dt = %.2f), shape = %s" %  ( str( data_specs ), t1, (t1-t0), str(variable.shape) ) )
        return variable, data_specs

    def loadFileFromCollection( self, collection, id=None ):
        t0 = time.time()
        wpsLog.debug( "loadFileFromCollection: '%s' '%s' " % ( collection, id ) )
        rv = self.collectionManager.getFile( collection, id )
        t1 = time.time()
        wpsLog.debug( "Done loadFileFromCollection, cdms2.open completed in %.2f sec" % (t1-t0) )
        return rv

if __name__ == "__main__":
    wpsLog.addHandler( logging.StreamHandler(sys.stdout) )
    wpsLog.setLevel(logging.DEBUG)

    dfile = 'http://dataserver.nccs.nasa.gov/thredds/dodsC/bypass/CREATE-IP/MERRA/mon/atmos/hur.ncml'
    slice_args = {'lev': (100000.0, 100000.0, 'cob')}
    dataset = f=cdms2.open(dfile)
    print dataset.id, dataset._status_
    dset = dataset( "hur", **slice_args )
    print str(dset.shape)

    slice_args1 = {"longitude": (-10.0, -10.0, 'cob'), "latitude": (10.0, 10.0, 'cob'), 'lev': (100000.0, 100000.0, 'cob')}
    dset1 = dataset( "hur", **slice_args1 )
    print str(dset1.shape)

    # from modules.configuration import MERRA_TEST_VARIABLES
    # test_point = [ -137.0, 35.0, 85000 ]
    # test_time = '2010-01-16T12:00:00'
    # collection = MERRA_TEST_VARIABLES["collection"]
    # varname = MERRA_TEST_VARIABLES["vars"][0]
    #
    # cache_type = CachedVariable.getCacheType( True, None )
    # region =  { "level":test_point[2] }
    # data = { 'name': varname, 'collection': collection }
    # variable, result_obj = dataManager.loadVariable( data, region, cache_type )
    # cached_data = dataManager.cacheManager.getVariable( result_obj['cache_id'], region )
    # cached_domain = cached_data[1]
    # cache_id = cached_domain.persist()
    # cached_domain.load_persisted_data()



