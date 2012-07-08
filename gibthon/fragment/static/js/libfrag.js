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

    //Fetch all available fragments, returned in a list
	this.getAll = function(_suc)
	{
       AJAX.post({
           url: '/fragment/api/listAll/',
           success: function(data)
           {
               var frags = new Array();
               for(f in data)
                {
                    frags.push(new Fragment(data[f].id, data[f].name, 
                                            data[f].desc, data[f].length));
                }
                _suc(frags);
           },
       });
	}
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

// ============================================ Widgets

/*
 *  jFragment: draggable fragment
 *
 * */
$.widget("ui.jFragment", $.ui.draggable, {
    options: {
        fragment: null,
        color: 'red',
    },
    _create: function() {
        this.f = this.options.fragment;
        this.$el = $(this.element[0]).draggable(this.options);
        this.$el.html("<p class='jf-name'>"+this.f.getName()+"</p>");
        this.$el.addClass('jFragment');
        this.$el.css({'background-color':this.options.color,
                     'border-color':this.options.color,});
    }
});  

/*
 *
 * jFragmentSelector: Easily select fragment from the available using drag and
 * drop
 *
 * */
$.widget("ui.jFragmentSelector", {
    options: {
        droptarget: null,
        containment: 'parent',
        dragEnabled: true,
    },
    _create: function() {
        var self = this;
        var $e = this.$el = $(this.element[0]);
        this.$fragView = $e.find('.JFS_fragView');
        this.$filterHolder = $e.find('.JFS_filterHolder');
        this.$filterHint = this.$filterHolder.find('p');
        this.$filterInput = this.$filterHolder.find('input');
        this.$numItems = $e.find('#JFS_num_items');
        //if we click the hint, we select the input
        this.$filterHint.on('click', function() {
            self.$filterInput.focus();
        });
        this.$filterInput.on('focus', function() {
            self._onInputFocus();
        });
        this.$filterInput.on('blur', function() {
            self._onInputBlur();
        });
        this.$filterInput.hover(function() {
            self.$filterHolder.addClass('ui-state-hover');   
        }, function() {
            self.$filterHolder.removeClass('ui-state-hover');   
        });
        console.log('jFragmentSelect _created');

        //fetch the fragments
        libFrag.getAll(function(frags){
            //remove the loading screen
            self.$fragView.empty();
            //add in the fragments one by one
            for(f in frags)
            {
                $('<div/>').addClass('JFS_fragHolder')
                .append( $('<div/>').jFragment({
                    fragment: frags[f], 
                    helper:'clone'
                }))
                .appendTo(self.$fragView);
            }
            self.$numItems.text(frags.length);
        });


    },
    //set height to fill the parent container
    height: function(h){
        console.log('jFragmentSelect height('+h+')');
        this.$fragView.height( h - 
                              (this.$el.height() - this.$fragView.height()));        
    },
    _onInputFocus: function(){
        this.$filterHint.hide();
        this.$filterHolder.addClass('ui-state-focus');
    },
    _onInputBlur: function(){
        if(!this.$filterInput.val())
            this.$filterHint.show();
        this.$filterHolder.removeClass('ui-state-focus');
    },

});
      
