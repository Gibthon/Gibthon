from BeautifulSoup import BeautifulSoup
import datetime
import urllib2
import re

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

'''
Stem of the URL for a BioBrick page
e.g. the xml info for the part "BBa_I13520" is found at:
"http://partsregistry.org/cgi/xml/part.cgi?part=BBa_I13520"
'''
#the url from which to fetch parts
partsURL = "http://partsregistry.org/cgi/xml/part.cgi?part=%s"

#regexes for all possible parts - add more if you discover them!
name_forms = [	"^bba_[a-z][0-9]{0,15}$",
				"^psb[0-9a-z]{0,15}$"]

def isValidPart(name):
	"""
	Check if name is a valid part name
	"""
	for n in name_forms:
		if re.match(n, name.lower()) is not None:
			return True
	return False	

class Part:
	def __init__(self, name):
		#see if the part name is sane
		if not isValidPart(name):
			raise ValueError("Invalid part name '%s'" % name)
		#try and retrieve BioBrick data
		try:
			print("Fetching data from:\t %s" % (partsURL % name) )
			f = urllib2.urlopen(partsURL % name)
			soup = BeautifulSoup(f)
			if soup.find("error"):
				raise ValueError("Part '%s' does not exist" % name)
			else:
				print "Successfully fetched part %s!" % name
		except urllib2.URLError:
			raise IOError("Could not retrieve data.")
		def getTag(asoup, tag, default=u""):
			if asoup.find(tag).contents:
				return asoup.find(tag).string.strip()
			else:
				return default
		#name, description, type, status, etc
		self.name = getTag(soup, "part_name")
		self.short_name = getTag(soup, "part_short_name")
		self.description = getTag(soup, "part_short_desc")
		self.status = getTag(soup, "part_status")
		self.results = getTag(soup, "part_results")
		self.nickname = getTag(soup, "part_nickname")
		self.rating = getTag(soup, "part_rating")
		self.url = getTag(soup, "part_url")
		self.author = getTag(soup, "part_author")
		self.best_quality = getTag(soup, "best_quality")
		#type and categories
		self.part_type = getTag(soup, "part_type")
		self.categories = [x.string.strip() for x in soup.findAll('category')]
		#date added
		adate = re.findall("[0-9]+", soup.find('part_entered').string)
		self.date_entered = datetime.date(int(adate[0]), int(adate[1]), int(adate[2]))
		#deep subparts
		self.deep_subparts = [x.string.strip() for x in soup.find('deep_subparts').findAll('part_name')]
		#specified subparts
		self.specified_subparts = [x.string.strip() for x in soup.find('specified_subparts').findAll('part_name')]
		#features
		self.features = list()
		for element in soup.findAll('feature'):
			direction = getTag(element, 'direction', u"forward").lower()
			if direction == 'forward':
				strand = 1
			elif direction == 'reverse':
				strand = -1
			else:
				strand = 0
			feature = { "name" : getTag(element, 'title', u"NoName"),
						"type" : getTag(element, 'type', u"NoType"),
						"strand": strand,
						"startpos" : int(element.find('startpos').string.strip()),
						"endpos" : int(element.find('endpos').string.strip())}
			self.features.append(feature)
		#sequence
		self.sequence = re.sub(r'\s', '', soup.find('seq_data').string)
		#twins
		self.twins = list()
		for element in soup.findAll('twin'):
			self.twins.append(element.string.strip())
		#scars
		#the names 'RFC[10]' and 'RFC[10.1]' are Randy's, I don't know what they mean
		for match in re.finditer("tactagag", self.sequence):
			ends = [x['endpos'] for x in self.features]
			starts = [x['startpos'] for x in self.features]
			if match.start() in ends and (match.end() + 1) in starts:
				self.features.append( { "name" : "RFC[10]",
										"type" : "scar",
										"strand": 1,
										"startpos" : match.start() + 1,
										"endpos" : match.end() } )
		for match in re.finditer("tactag", self.sequence):
			ends = [x['endpos'] for x in self.features if x['type'].lower() == u'rbs']
			starts = [x['startpos'] for x in self.features]
			if match.start() in ends and (match.end() + 1) in starts:
				self.features.append( { "name" : "RFC[10.1]",
										"type" : "scar",
										"strand": 1,
										"startpos" : match.start() + 1,
										"endpos" : match.end() } )
		#close the filehandle
		f.close()
		
	def to_seq_record(self):
		"""Return a SeqRecord object containing all the data from the Part"""
		#create the anotations in a pythonic manner
		exempt = ['name', 'description', 'features', 'sequence'] #things which aren't annotations
		annotations = { }
		for key, value in self.__dict__.iteritems():
			if key.lower() not in exempt:
				annotations[key] = value
		
		#create the features
		features = []
		for feat in self.features:
			features.append( SeqFeature( 
				location = FeatureLocation(feat['startpos'], feat['endpos']),
				type = feat['type'],
				strand = feat['strand'],
				id = feat['name']))
		
		return SeqRecord(	self.sequence, 
							id=self.name,
							name=self.name,
							description=self.description,
							features=features,
							annotations=annotations)
		
