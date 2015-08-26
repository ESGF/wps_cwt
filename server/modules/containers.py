import json
from utilities import  *

class JSONObject:

    def __init__( self, spec={} ):
        if isinstance( spec, JSONObject ): spec = spec.spec
        self.spec = convert_json_str( spec )
        self.items = {}
        assert isinstance( self.spec, dict ), "Error, unrecognized JSONObject spec: %s " % str( spec)
        self.process_spec( )

    def process_spec( self ):
        self.items = dict( self.spec )

    def __str__( self ):
        return dump_json_str(self.items)

    def get(self, key, default_val=None ):
        return self.items.get( key, default_val )

    def items(self):
        return self.items.items()

    def iteritems(self):
        return self.items.iteritems()

    def __getitem__(self, item):
        return self.items.get( item, None )

    def __setitem__(self, key, value):
        self.items[ key ] = value

    def getItem( self, feature_name ):
        return self.items.get( feature_name, None )

class JSONObjectContainer:

    def __init__( self, spec=None ):
        if isinstance( spec, JSONObjectContainer ): spec = spec.spec
        self._objects = []
        self.process_spec( spec )

    def newObject( self, spec ):
        return JSONObject(spec)

    def process_spec( self, spec ):
        if spec:
            spec = convert_json_str( spec )
            if isinstance( spec, list ):
                for object_spec in spec:
                    self._objects.append( self.newObject( object_spec) )
            elif isinstance( spec, dict ):
                self._objects.append( self.newObject( spec) )
            else:
                raise Exception( "Unrecognized JSONObject spec: " + str(spec) )

    @property
    def value(self):
        return self._objects[0] if len(self._objects) else None

    @property
    def values(self):
        return self._objects