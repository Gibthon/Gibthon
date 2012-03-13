/*
 * libfrag.js -- Javascript API for the fragment library
 * --jQuery must already have been imported
 * * * * * * * * * * "Public" Functions - 
 * 
 * list_fragments( func success([metadata]) )
 * 		get the metadata of all the fragments available to the user
 * 
 * get_meta(fid, func success(metadata))
 * 		get the metadata of a particular fragment
 * 
 * set_meta(fid, metadata, func success())
 * 		set the metadata of a particular fragment
 * 
 * get_feats(fid, func success(data))
 * 		get a particular fragment's features, data.feats, its length is data.length and data.alpha is dictionary containing
 * 		the alphabet
 * 
 * set_feats(fid, features, func success())
 * 		set a particular fragment's features
 * 
 * get_length(fid)
 * 		return the length of the sequence
 * 
 * get_sequence(fid, start, length, func success(sequence))
 * 		returns a string containing a section of the sequence, starting at start and of length length
 * 		if the sequence is shorter than start + length, the return sequence will be shorter too
 * 		if length is zero or negative, the whole sequence is fetched
 * 
 * * * * * * * * * * * Datatypes
 * 
 * metadata:
 * 		.id	: Fragment ID
 * 		.name	: Fragment Name
 *		.desc	: Description of fragment
 *		.refs	: A list of references, see below
 *		.annots	: Annotations, a list of [key, value] 'tuples'
 *		.origin	: Where the fragment came from, is not updated with set
 *		.length	: how long the sequence is, also not updated with set
 * 
 * reference:
 * 		.title		: Reference title
 *		.authors	: Authors
 *		.journal	: Journal
 *		.medline_id	: Blank if there isn't one
 *		.pubmed_id	: Blank if there isn't one
 * 
 * feature:
 * 		start		: Feature start position
 * 		end			: Feature end position
 * 		strand		: 1 or -1
 * 		type		: String
 * 		id			: 
 * 		qualifiers	: A list of qualifiers
 * 
 * qualifier:
 * 		name	: Qualifier Name
 * 		value	: value
 * 
 * sequence:
 * 		.seq			: A simple string containing the whole sequence
 * 		.complement()	: Return the sequence's complement
 * 		.rcomplement()	: Return the sequence's reverse complement
 * 		.reverse()		: Return the sequence in reverse
 * 
 * * * * * * * * * * * Other Functions
 * complement(seq, alphabet)
 * 		Return the complement of a sequence as a string
 * 
 * rcomplement(seq, alphabet) 
 * 		return the reverse complement of a sequence as a string
 * 
 * */

var PRE = '/fragment/api/'

var alphabet = {
	'A': 'T',
	'B': 'V',
	'C': 'G',
	'D': 'H',
	'G': 'C',
	'H': 'D',
	'K': 'M',
	'M': 'K',
	'N': 'N',
	'R': 'Y',
	'S': 'S',
	'T': 'A',
	'V': 'B',
	'W': 'W',
	'Y': 'R',
	'a': 't',
	'b': 'v',
	'c': 'g',
	'd': 'h',
	'g': 'c',
	'h': 'd',
	'k': 'm',
	'm': 'k',
	'n': 'n',
	'r': 'y',
	's': 's',
	't': 'a',
	'v': 'b',
	'w': 'w',
	'y': 'r',
}

// =========================== Other Functions
var complement = function(s, a)
{
	for(var i = 0; i < s.length; i = i + 1)
		s[i] = a[s[i]];
	return s;
}

var rcomplement = function(s, a)
{
	return complement(s.split("").reverse().join(""));
}

// ============================Constructors


/*
 * Represents a fragment and associated metadata 
 * */
function Fragment(d) // d = json from server
{
	var self = this;
	this.name = d.name;
	this.fid = d.fid;
	this.id = d.fid;
	this.desc = d.desc;
	this.refs = new Array();
	for(var i in d.refs)
	{
		this.refs.push(new Reference(d.refs[i]));
	}
	this.annots = new Array();
	for(var i in d.annots)
	{
		this.annots.push([i, d.annots[i]]);
	}
	this.feats = new Array();
	for(var i in d.feats)
	{
		this.feats.push(new Feature(d.feats[i]));
	}
	this.origin = d.origin;
	this.length = d.length;
	
	this.seq = null;
	
	this.getFeatById = function(id)
	{
		for(var i in self.feats)
		{
			if(self.feats[i].id == id)
				return self.feats[i];
		}
		return null;
	}
		
}

function Reference(d)
{
	this.title = d.title;
	this.authors = d.authors;
	this.journal = d.journal;
	this.medline_id = d.medline_id;
	this.pubmed_id = d.pubmed_id;
}

function Feature(d)
{
	this.start = d.start;
	this.end = d.end;
	this.strand = d.strand;
	this.type = d.type;
	this.id = d.id;
	this.strand = d.strand;
	this.qualifiers = new Array();
	for(var i in d.qualifiers)
	{
		this.qualifiers.push(new Qualifier(d.qualifiers[i]));
	}
}

function Qualifier(d)
{
	this.name = d.name;
	this.value = d.value;
}

/*
 * 		Constructors depreciated from here on in
 * */
function Metadata(fid, name, desc, refs, annots, origin, length)
{
	this.fid = fid;
	this.name = name;
	this.desc = desc;
	this.refs = refs;
	this.annots = annots;
	this.origin = origin;
	this.length = length;
}

function Sequence(s)
{
	self = this;
	this.seq = s;
	this.complement = function() {return complement(self.seq, alphabet);}
	this.rcomplement = function() {return rcomplement(self.seq, alphabet);}
	this.reverse = function() {return self.seq.split("").reverse().join("");}
}
// ========================================================= AJAX things
function handle_error(action, textStatus, errorThrown)
{
	var s = "Error while " + action + ": \n\tstatus: " + textStatus + "\n\tthrew: " + errorThrown;
	console.error(s);
}

function make_request(url, data, desc, cb, error_cb)
{
	$.ajax({
		url: url,
		dataType: 'json',
		data: data,
		type:'POST',
		error: function(jqXHR, textStatus, errorThrown)
		{
			if(error_cb == undefined)
				handle_error(desc, "Ajax request failed, status: '" + textStatus + "'", errorThrown);
		},
		success: function(data) 
		{
			if(data[0] != 0)
			{
				if(error_cb == undefined)
					handle_error(desc, "Error at Server", data[1]);
				return;
			}
			cb(data[1]);
		},
	});
}

// ======================================== library functions

var list_fragments = function(success)
{
	make_request(PRE + 'listAll/', undefined, "getting fragment list", success);
}

var get_meta = function(fid, success)
{
	make_request(PRE + fid + '/getMeta/', undefined, "getting metadata for '" + fid +"'", success);
}

var get_multi_meta = function(fids, success)
{
	var d = JSON.stringify({'fids[]': fids,});
	make_request(PRE + '/getMeta/', d, "getting metadata for '" + fids +"'", success);
}

var set_meta = function(fid, metadata, success)
{
	make_request(PRE + fid + '/setMeta/', JSON.stringify(metadata), "saving metadata for '" + fid +"'", success);
}

var get_feats = function(fid, success)
{
	make_request(PRE + fid + '/getFeats/', undefined, "getting features for '" + fid +"'", success);
}

var set_feats = function(fid, features, success)
{
	make_request(PRE + fid + '/setFeats/', JSON.stringify(features), "saving features for '" + fid +"'", success);
}

var get_length = function(fid, success)
{
	make_request(PRE + fid + '/getLength/', undefined, "getting length for '" + fid +"'", success);
}

var get_sequence = function(fid, start, length, success)
{
	make_request(PRE + fid + '/getSeq/', undefined, "getting sequence for '" + fid +"' offset:" + start + " length:" + length, success);
}
