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
 * get_feats(fid, func success([features]))
 * 		get a particular fragment's features
 * 
 * set_feats(fid, features, func success())
 * 		set a particular fragment's features
 * 
 * get_sequence(fid, start, length, func success(sequence))
 * 		returns a string containing a section of the sequence, starting at start and of length length
 * 		if the sequence is shorter than start + length, the return sequence will be shorter too
 * 		if length is zero or negative, the whole sequence is fetched
 * 
 * * * * * * * * * * * Datatypes
 * 
 * metadata:
 * 		.fid	: Fragment ID
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
 *		.alpha			: 'Dictionary' containing all possible nucleotides as keys, and all complements as values
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

// =========================== Other Functions
var complement(s, a)
{
	for(var i = 0; i < s.length; i = i + 1)
		s[i] = a[s[i]];
	return s;
}

var rcomplement(s, a)
{
	return complement(s.split("").reverse().join(""));
}

// ============================Constructors
var Metadata(fid, name, desc, refs, annots, origin, length)
{
	this.fid = fid;
	this.name = name;
	this.desc = desc;
	this.refs = refs;
	this.annots = annots;
	this.origin = origin;
	this.length = length;
}

var Reference(title, authors, journal, medline_id, pubmed_id)
{
	this.title = title;
	this.authors = authors;
	this.journal = journal;
	this.medline_id = medline_id;
	this.pubmed_id = pubmed_id;
}

var Feature(start, end, strand, type, id, qualifiers)
{
	this.start = start;
	this.end = end;
	this.strand = strand;
	this.type = type;
	this.id = id;
	this.qualifiers = qualifiers;
	if(this.start < end) this.strand = 1;
	if(this.start > end) this.strand = -1;
}
 
var Qualifier(name, value)
{
	this.name = name;
	this.value = value;
}

var Alphabet(s, a)
{
	self = this;
	this.seq = s;
	this.a = a;
	this.complement = function() {return complement(self.seq, self.a);}
	this.rcomplement = function() {return rcomplement(self.seq, self.a);}
	this.reverse = function() {return self.seq.split("").reverse().join("");}
}
// ========================================================= AJAX things
var handle_error(action, msg)
{
	var s = "Error while " + action + ": " + msg;
	//alert(s);
	console.error(s);
}

var make_request(u, d, act, cb)
{
	$.ajax({
		url: u,
		dataType: 'json',
		data: d,
		type:'POST',
		error: function(jqXHR, textStatus, errorThrown)
		{
			handle_error(act, "Ajax request failed, status: '" + textStatus + "', error: '" + errorThrown + "'"
		}
		success: function(data) 
		{
			if(data[1] != 0)
			{
				handle_error(act, data[0]);
				return;
			}
			cb(data[0]);
		}
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

var set_meta = function(fid, metadata, success)
{
	make_request(PRE + fid + '/setMeta/', {'meta':metadata,}, "saving metadata for '" + fid +"'", success);
}

var get_feats = function(fid, success)
{
	make_request(PRE + fid + '/getFeats/', undefined, "getting features for '" + fid +"'", success);
}

var set_feats = function(fid, features, success)
{
	make_request(PRE + fid + '/setFeats/', {'features':features,}, "saving features for '" + fid +"'", success);
}

var get_sequence = function(fid, start, length, success)
{
	make_request(PRE + fid + '/getSeq/', undefined, "getting sequence for '" + fid +"' offset:" + start + " length:" + length, success);
}
