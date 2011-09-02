/*
formExtender:
* expands a single form field into multiple form fields
* 
* expected initialHTML
<div id='extender'> // the form extender
	<div id='extender-copy' style='display:none'> 
		everything in here is copied
		<button class='extender-remove'></button>
	</div>
	<div id='extender-target'></div>
	
	<button class='extender-add'></button>
</div>
*/

(function( $, undefined ) {

$.widget("ui.formExtender", {
	options: {
		addInitial: true, //whether to add an initial element
	},
	_init: function() {
		var self = this;
		this.$el = $(this.element[0]);
		
	//intialise the elements	
		this.$copy = this.$el.find('#extender-copy').detach(); //remove from the dom - it'll only cause trouble
		this.$target = this.$el.find('#extender-target');
		
		this.$el.find('button.extender-add').button({
			label: 'More',
			icons: {primary: 'ui-icon-plusthick',},
		}).click( function(event) {
			self._add(event);
		});
	//intialise the data
		this.num_elements = 0;
		
		if(this.options.addInitial)
		{
			this._add();
		}
	},
	numElements: function() {return num_elements},
	_add: function(event){
		var self=this;
		var $clone = this.$copy.clone();
		
		$clone.find('button.extender-remove').button({
			text: false,
			icons: {primary: 'ui-icon-trash',},
		}).click( function(event) {
			console.log('remove');
			self._remove(event);
		});
		
		$clone.appendTo(this.$target)
			.slideDown('fast');
		this.num_elements = this.num_elements + 1;
	},
	_remove: function(event){
		if(this.num_elements <= 1)
			return;
		$(event.target).closest('#extender-copy').slideUp('fast', function() {
			$(this).remove();
		});
		this.num_elements = this.num_elements - 1;
	},
});

})(jQuery);

/*
 * magicForm
 * - forms which appear when needed and submit&update via AJAX.
 * Expected HTML:
<form action,method,etc>
	<button class='magic-button'> //the magic button!
	<div class='magic-item'>
		<p class='magic-text'></p>
		<input(s) name='name' class='magick-input'/>
		<p class='magic-error'></p>
	</div>
</form>
 */
(function( $, undefined ) {

$.widget("ui.magicForm", {
	options: {
		save: function(data) {}, //called on save
		edit: function() {},
	},
	_init: function() {
		var self = this;
		this.$form = $(this.element[0]);
		this.action = this.$form.attr('action');
		this.method = this.$form.attr('method');
		if((this.method == '') || (this.method == undefined))
			this.method = 'GET'
		
		this.edit_opts = {'label': 'Edit', 'icons': {'primary': 'ui-icon-pencil',},}
		this.save_opts = {'label': 'Save', 'icons': {'primary': 'ui-icon-disk',},}
		this.$button = this.$form.find('button.magic-button')
			.button(this.edit_opts)
			.click( function() {self._edit();});
	},
	method: function(method)
	{
		if(method == undefined)
			return this.method;
		this.method = method;
	},
	action: function(action)
	{
		if(action == undefined)
			return this.action;
		this.action = action;
	},
	_edit: function()
	{
		var self = this;
		this.$form.find('.magic-item').each( function() {
			var $item = $(this);
			var $input = $item.find('.magic-input');
			var $text = $item.find('.magic-text');
			//set the value
			$input.val($text.text());
			$text.fadeOut('fast', function() {$input.fadeIn('fast');});
		});
		this.$button.button('option', this.save_opts)
			.unbind('click')
			.click( function() {self._save();});
		this._trigger('edit');
	},
	_save: function()
	{
		var self = this;
		$.ajax({
			url: this.action,
			dataType: 'json',
			data: this.$form.serialize(),
			type: this.method,
			success: function(data) {self._data(data);},
		});
	},
	_data: function(data)
	{
		var self = this;
		if(data[0] != 0) //if error
		{
			console.log('Error!');
			$inputs = this.$form.find('.magic-input');
			for(i in data[1].errors)
			{
				console.log('errors.' + i + ' = ' + data[1].errors[i]);
				//find the input(s) to which the error relates
				$inputs.each( function() {
					var $t = $(this);
					if($t.attr('name') == i)
					{
					console.log('input name: ' + $t.attr('name') + ' i: ' + i);
						//show the error
						var $item = $t.closest('.magic-item').find('.magic-error');
						console.log('  item.html(): ' + $item.html());
						console.log('  errors[i]: ' + data[1].errors[i]);
						$item.text('' + data[1].errors[i]);
						$item.slideDown('slow');
					}
				});
			}
		}
		else
		{	//hide the errors
			this.$form.find('.magic-error').each(function() {$(this).slideUp('slow');});
			this.$button.button('option', this.edit_opts)
				.unbind('click')
				.click( function() {self._edit();});
			this.$form.find('.magic-item').each( function() {
				var $item = $(this);
				var $input = $item.find('.magic-input');
				var $text = $item.find('.magic-text');
				//show the text
				$input.fadeOut('fast', function() {
					var val = 'None';
					if(data[1].fields != undefined)
						val = data[1].fields[$input.attr('name')];
					$text.text(val).fadeIn('fast');
				});
			});
		}
	},
});

})(jQuery);
