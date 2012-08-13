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
				_suc(new Fragment(data));
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
                    frags.push(new Fragment(data[f]));
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
function Fragment(data)
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
     * cropsettings = {
     *  start: start position,
     *  end: end position,
     *  f_internal: keep internal features?,
     *  f_all: keep all features
     *  result: in ['new', 'overwrite']
     *  new_name: defined if result == new
     *  new_desc: defined if result == new
     *  };
     *
     *  */
    this.crop = function(cropsettings, cb)
    {
        AJAX.post({
            url: '/fragment/api/'+ id + '/crop/',
            data: cropsettings,
            success: function(data){
                if($.isFunction(cb)) cb(data);
            },
            error: function() {
                console.error('Error cropping fragment');
            },
        });
    };

    this.toString = function()
    {
        return '[Fragment name="' + name + '"]';
    }

	/*
	 *
	 * Private data
	 *
	 */

	var id = data.id;
	var name = data.name;
	var desc = data.desc;
	var length = data.length;

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
        draggable: true,
    },
    _create: function() {
        this.f = this.options.fragment;
        this.$el = $(this.element[0])
        if(this.options.draggable)
            this.$el.draggable(this.options);
        this.$el.html("<p class='jf-name'>"+this.f.getName()+"</p>");
        this.$el.addClass('jFragment');
        this.$el.css({'background-color':this.options.color,
                     'border-color':this.options.color,});
    },
    getFragment: function()
    {
        return this.f;
    },
    getName: function()
    {
        return this.f.getName();
    },
    getColor: function()
    {
        return this.$el.css('background-color');
    },
    option: function(name, value)
    {
        if(value == undefined)
            return this.options[name];
        console.log('setting '+name);
        this.options[name] = value;
    },
});  

/*
 *
 * jFragmentSelector: Easily select fragment from the available ones 
 * using drag and drop
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

        this.filter_timeout = null;
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
        this.$filterInput.keypress( function(){
            if(self.filter_timeout != null){
                clearTimeout(self.filter_timeout);
                self.filter_timeout = null;
            }
            self.filter_timeout = setTimeout(function(){
                self._filter();
            }, 350);
        });
        this.$filterInput.hover(function() {
            self.$filterHolder.addClass('ui-state-hover');   
        }, function() {
            self.$filterHolder.removeClass('ui-state-hover');   
        });

        this.frags = new Array();
        //fetch the fragments
        libFrag.getAll(function(frags){
            //remove the loading screen
            self.$fragView.empty();
            //add in the fragments one by one
            for(f in frags)
            {
                self.frags.push(frags[f]);
                $('<div/>').addClass('JFS_fragHolder')
                .append( $('<div/>').jFragment({
                    fragment: frags[f], 
                    helper: function(){
                        return $('<div/>').jFragment({
                            draggable:false,
                            fragment: $(this).jFragment('getFragment'),
                        }).appendTo(self.$el);
                    },
                    containment: self.options.containment,
                    zIndex:200,
                }))
                .data('f', f)
                .appendTo(self.$fragView);
            }
            self.$numItems.text(frags.length);
        });


    },
    //set height to fill the parent container
    height: function(h){
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
    _filter: function(){
        var s = this.$filterInput.val();
        var self = this;
        var matches = 0;
        if(s){
            this.$el.find('.JFS_fragHolder').each(function(index, el){
                //hide the element if it doesn't match s
                if( self.frags[index].getName().toLowerCase()
                   .indexOf(s.toLowerCase()) < 0)
                {
                    $(this).hide();
                }
                else
                {
                    $(this).show();
                    matches = matches + 1;
                }
            });
        }
        else{
            this.$el.find('.JFS_fragHolder').each(function(index, el){
                //show all if s is empty
                $(this).show();
                matches = matches + 1;
            });
        }
        this.$numItems.text(matches);
    },
    option: function(name, value)
    {
        if(value == undefined)
            return this.options[name];
        console.log('setting '+name);
        this.options[name] = value;
    },
});
      
