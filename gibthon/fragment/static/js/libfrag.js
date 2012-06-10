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
	this.getByID = function(id, success)
	{
		AJAX.post({
			url: '/api/'+id+'/', 
			success: function(data)
			{
				return new Fragment(data.id, data.name, data.desc, data.length);
			},
		});
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
	var getID = function() {return id;};
	var getName = function() {return name;};
	var getDesc = function() {return desc;};
	var getLength = function() {return length;};

	//if the sequence has not already been fetched, it will be streamed, otherwise
	//complete_fn will be called immediately
	var getSequence = function( update_fn, complete_fn)
	{
		if(sequence!=null)
		{
			complete_fn(sequence);
			return;
		}

		//stream the sequence from the server
		AJAX.stream({
			url: 'api/' + id + '/getSeq/',
			success: function(data)
			{
				sequence=data;
				complete_fn(sequence);
			},
		},
		update_fn);

	}

	//if the metadata has not already been fetched, it will be fetched otherwise
	//complete_fn will be called immediately
	var getMeta = function(complete_fn)
	{
		if(metadata!=null)
		{
			complete_fn(metadata);
			return;
		}

		//fetch the metadata from the server
		AJAX.post({
			url: 'api/' + id + '/getMeta/',
			success: function(ret){
				meta = ret;
				complete_fn(ret);
			},
			error: function() {},
		});
	}

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

