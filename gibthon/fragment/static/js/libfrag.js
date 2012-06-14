/*
 * Fragment API
 *  --fetch fragment from the server
 *
 */

var libFrag = new function()
{
	/*
	 * Public functions
	 *
	 */

	//fetch a fragment by its ID
	this.getByID = function(id, _suc)
	{
		AJAX.post({
			url: '/fragment/api/'+id+'/', 
			success: function(data)
			{
				_suc(new Fragment(data.id, data.name, data.desc, data.length));
			},
		});
	};

	//IUAPC Unambiguous Alphabet / Complement table
	this.alphabet = {
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
		' ': ' ',
	}
}

/*
 *
 * An actual fragment object
 *
 */
function Fragment(_id, _name, _desc, _length)
{
	/*
	 * Public Accessors
	 *
	 */
	this.getID = function() {return id;};
	this.getName = function() {return name;};
  this.getDesc = function() {return desc;};
	this.getLength = function() {return length;};

	//if the sequence has not already been fetched, it will be streamed, otherwise
	//complete_fn will be called immediately
	this.getSequence = function( update_fn, complete_fn)
	{
		if(sequence!=null)
		{
			complete_fn(sequence);
			return;
		}

		//stream the sequence from the server
		AJAX.stream({
			url: '/fragment/api/' + id + '/getSeq/',
			success: function(data)
			{
				sequence = new String();
				for(i in data)
				{
					if(data[i] != ' ')
						sequence = sequence + data[i];
				}
				complete_fn(sequence);
			},
			type: 'GET',
			error: function(error)
			{
				console.error('Error fetching sequence: ' + error);
			}
		},
		function(data)
		{
			sequence = new String();
			for(i in data)
				{
					if(data[i] != ' ')
						sequence = sequence + data[i];
				}
			update_fn(sequence);
		} );
	};

	//if the metadata has not already been fetched, it will be fetched otherwise
	//complete_fn will be called immediately
	this.getMeta = function(complete_fn)
	{
		if(metadata!=null)
		{
			complete_fn(metadata);
			return;
		}

		//fetch the metadata from the server
		AJAX.post({
			url: '/fragment/api/' + id + '/getMeta/',
			success: function(ret){
				meta = ret;
				complete_fn(ret);
			},
			error: function() {},
		});
	};

	this.setMeta = function(new_meta, success_cb, fail_cb)
	{
		metadata = new_meta;
		AJAX.post({
			url: '/fragment/api/' + id + '/setMeta/',
			success: success_cb,
			error: fail_cb,
			data: metadata,
		});
	};	

	this.getFeats = function(success_fn)
	{
		AJAX.post({
			url: '/fragment/api/' + id + '/getFeats/',
			success: success_fn,
			error: function(err)
			{
				console.error('Error getting features for ' + name);
			},
		});
	};

	/*
	 *
	 * Private data
	 *
	 */

	var id = _id;
	var name = _name;
	var desc = _desc;
	var length = _length;

	//Data which might get filled in by subsequent AJAX calls
	var sequence = null;
	var metadata = null;
}
