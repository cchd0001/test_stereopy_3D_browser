import re,os,time
from http.server import BaseHTTPRequestHandler, HTTPServer
from typing import Optional, Union
import anndata as ad
from anndata import AnnData

class Stereo3DWebCache:
    """
    Analyse the 3D SRT data and provide detailed json data for the data browser.
    """
    def __init__(self,
                 adata: AnnData,
                 meshes: {},
                 cluster_label:str = 'Annotation',
                 spatial_label:str = 'spatial_rigid',
                 geneset = None):
        self._data = adata
        self._annokey = cluster_label
        self._spatkey = spatial_label
        self.InitMeshes(meshes)
        
    def get_summary(self):
        """
        return the summary.json
        """
        return ""

    def get_gene(self,genename):
        """
        return the Gene/xxxgene.json
        """
        return ""
    
    def get_genename(self):
        """
        return the gene.json
        """
        return ""
    
    def get_mesh(self,meshname):
        """
        return the meshes.json
        """
        return ""
    
    def get_anno(self,annoname):
        """
        return the Anno/xxxanno.json
        """
        return ""
    

class StoppableHTTPServer(HTTPServer):
    """
    The http server that stop when not_forever is called.
    """
    def serve_forever(self):
        self.stopped = False
        while not self.stopped:
            self.handle_request()
            time.sleep(0.100)
            
    def not_forever(self):
        print('Server terminate ...',flush=True)
        self.stopped = True
        self.server_close()
        
class ServerDataCache:
    def __init__(self):
        self._data_hook = None
        self._server = None
    
    @property
    def data_hook(self):
        return self._data_hook
    @property
    def server(self):
        return self._server
    
    @server.setter
    def server(self,http):
        self._server = http
    
    @data_hook.setter
    def data_hook(self, data_hook):
        self._data_hook = data_hook
        

ServerInstance = ServerDataCache()

class DynamicRequstHander(BaseHTTPRequestHandler):
    """
    The request hander that return static browser files or detailed data jsons.
    """

    def stop_server(self):
        ServerInstance.server.not_forever()
        self.send_response(404)
        self.send_header('Content-type', 'text/html')
        self.end_headers()
        self.wfile.write(b"Server shotdown now!")
        
    def ret_static_files(self, the_relate_path, file_type):
        """
        return all static files like the browser codes and images
        """
        data_dir = "C:/Users/guolidong/Desktop/test"
        visit_path = f"{data_dir}{the_relate_path}"
        try:
            self.send_response(200)
            self.send_header('Content-type', file_type)
            self.end_headers()
            f = open(visit_path, 'rb')
            self.wfile.write(f.read())
            f.close()
        except:
            self.ret_404()
        
    def ret_404(self):
        """
        return 404
        """
        self.send_response(404)
        self.send_header('Content-type', 'text/html')
        self.end_headers()
        self.wfile.write(b"404 Not Found")
    
    def do_GET(self):
        print(f'{self.path}',flush=True)
        #if len(self.path.split('?')) > 1:
        #    self.args = self.path.split('?')[1]
        self.path = self.path.split('?')[0]
        # set index.html as root by default
        if self.path in ['','//','/','/index.html']:
            self.ret_static_files("/index.html", 'text/html')  
        elif self.path == '/endnow':
            self.stop_server()
        elif self.path == '/test.json':  #handle json requst in the root path
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            ret_json = '{"test01":1.1, "test02":[1.1,3,2]}'
            self.wfile.write(bytes(ret_json,'UTF-8'))
        else:
            # extract xxx from api/xxx using regex
            match_html = re.search('(.*).html$', self.path)
            match_js = re.search('(.*).js$', self.path)
            match_ttf = re.search('(.*).ttf$', self.path)
            match_API_json = re.search('/api/(.*).json', self.path)
            if match_js:
                self.ret_static_files(self.path , 'application/javascript')
            elif match_html:
                self.ret_static_files(self.path , 'text/html')
            elif match_ttf:
                self.ret_static_files(self.path ,'application/x-font-ttf')
            elif match_API_json:
                #print('----------')
                xxx = match_API_json.group(1)
                # call function ABC with xxx as parameter
                self.send_response(200)
                self.send_header('Content-type', 'application/json')
                self.end_headers()
                ret_json = '{"test01":"{' +xxx+'}\", \"test02\":[1.1,3,2,444444]}'
                self.wfile.write(bytes(ret_json,'UTF-8'))
                #print('----------')
            else: # may be this is an static local file ? let try it
                self.ret_404()


                
def launch(datas,
           meshes:{},
           port:int = 7654,
           cluster_label:str = 'Annotation',
           spatial_label:str = 'spatial_rigid',
           geneset = None):
    """
    Launch a data browser server based on input data
    
    :param datas: an AnnData object or a list of AnnData objects
    :param mesh: all meshes in dict like : {'heart': 'pathxxx/heart.obj', liver:'pathxxx/liver.obj'}
    :parma port: the port id
    :param cluster_label: the keyword in obs for cluster/annotation info
    :param spatial_label: the keyword in obsm for 3D spatial coordinate
    :param geneset: the specific geneset to display, show all genes in var if geneset is None
    
    :return:
    """
    #merge anndata if necessary
    if type(datas) == list:
        if len(datas) < 1:
            print('No data provided, return without any data browsing server...')
            return
        adata  = datas[0]
        if len(datas) > 1:
            for i in range (1,len(datas)):
                adata = adata.concatenate(h5datas[i])
    else:
        adata = datas
    #sanity check for parameters
    if not (cluster_label in adata.obs.columns and spatial_label in adata.obsm):
        print('invalid keyword provided, return without any data browsing server...')
        return
    for meshname in meshes:
        meshfile = meshes[meshname]
        if os.path.isfile(meshfile):
            print(f'invalid obj :{meshfile}, return without any data browsing server...')
            return
    #create core datacache
    datacache = Stereo3DWebCache(adata,meshes,cluster_label,spatial_label,geneset)
    ServerInstance.data_hook = datacache
    
    #create webserver
    server_address = ('', port)
    httpd = StoppableHTTPServer(server_address, DynamicRequstHander)
    ServerInstance.server = httpd
    #start endless waiting now...
    print(f'Starting server on http://127.0.0.1:{port}')
    print(f'To ternimate this server, please click the close button.\n or visit http://127.0.0.1:{port}/endnow to end this server.')
    httpd.serve_forever()
    
