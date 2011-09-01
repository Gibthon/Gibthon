import httplib
import json
import urllib

class MessagePasser():
	def __init__(self, _channel):
		self.channel = _channel
	
	def conn(self):
		return httplib.HTTPConnection("async-message-passer.appspot.com")
	
	def fetch(self):
		conn = self.conn()
		conn.request("GET", "/?channel_name=%s"%(self.channel))
		res = conn.getresponse()
		if res.status != 200:
			return False
		else:
			data = res.read()
			return json.loads(data)
			
	def clear(self):
		conn = self.conn()
		params = urllib.urlencode({'clear_channel':1})
		conn.request("POST", "/submit?channel_name=%s"%(self.channel), params)
		res = conn.getresponse()
		return res.status == 200
	
	def get_key(self):
		conn = self.conn()
		conn.request("GET", "/?channel_name=%s&response_format=plain_with_push"%(self.channel))
		res = conn.getresponse()
		data = res.read()
		print json.loads(data)