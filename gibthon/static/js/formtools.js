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
		unique: false, //whether to give each input an unique name
		addInitial: true, //whether to add an initial element
		disabled: false,
		beforeAdd: function(e, $el) {}, //called before adding an element to the DOM
	},
	_init: function() {
		var self = this;
		this.$el = $(this.element[0]);
		
	//intialise the elements	
		this.$copy = this.$el.find('#extender-copy').detach(); //remove from the dom - it'll only cause trouble
		this.$target = this.$el.find('.extender-target');
		
		this.$el.find('button.extender-add').button({
			label: 'Add',
			icons: {primary: 'ui-icon-plusthick',},
		}).click( function(event) {
			self._add(event);
		});
		
		this.$el.find('button.extender-remove').each( function() {
			$(this).button({
				text: false,
				icons: {primary: 'ui-icon-trash',},
			}).click( function(event) {
				console.log('remove');
				self._remove(event);
			});
		});
		
	//intialise the data
		this.num_elements = 0;
		if(this.options.unique)
		{
			$('.extender-item').each( function() {
				$(this).find('input').each( function () {
					$(this).attr('name', $(this).attr('name') + self.num_elements);
				});
				self.num_elements = self.num_elements + 1;
			});
		}
		
		if(this.options.addInitial)
		{
			this._add();
		}
		if(this.options.disabled) this.disable();
	},
	enable: function()
	{
		this._enable_disable('enable');
	},
	disable: function()
	{
		this._enable_disable('disable');
	},
	numElements: function() {return num_elements},
	_enable_disable: function(fn)
	{
		var self = this;
		this.$el.find('button.extender-remove').each( function() {
			$(this).button(fn);
		});
		this.$el.find('button.extender-add').button(fn);
	},
	_add: function(event){
		var self=this;
		var $clone = this.$copy.clone().addClass('extender-item');
		
		if(this.options.unique)
		{
			$clone.find('input').each( function() {
				$(this).attr('name', $(this).attr('name') + self.num_elements);
			});
			self.num_elements = self.num_elements + 1;
		}
		
		$clone.find('button.extender-remove').button({
			text: false,
			icons: {primary: 'ui-icon-trash',},
		}).click( function(event) {
			self._remove(event);
		});
		
		this._trigger('beforeAdd', event, $clone);
		
		$clone.appendTo(this.$target)
			.slideDown('fast');
		this._trigger('add', event, $clone);
	},
	_remove: function(event){
		var self = this;
		if(this.num_elements <= 1)
			return;
		$(event.target).closest('.extender-item').slideUp('fast', function() {
			$(this).remove();
			self._trigger('remove', event, undefined);
		});
	},
});

})(jQuery);

var httppat = /^http:\/\//i;

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
		save: function(ev, data) {}, //called on save
		cancel: function() {},//called on cancel
		edit: function() {},
		defaultInput: "<input class='magic-input' type='text'></input>",
		defaultError: "<p class='magic-error'></p>",
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
		this.$button_cancel = this.$form.find('button.magic-button-cancel')
			.button({label: 'Cancel', icons: {primary: 'ui-icon-cancel'}, disabled: true,})
			.click( function() {self._cancel();});
			
		this.$form.find('.magic-item').each( function() {
			var $text = $(this).find('.magic-text');
			//check whether there's an input
			if($(this).find('.magic-input').length == 0)
			{
				$(self.options.defaultInput).attr('name', $(this).attr('id')).insertAfter($text);
			}
			//check whether there's an error
			if($(this).find('.magic-error').length == 0)
			{
				$(this).append($(self.options.defaultError));
			}
		});
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
			$text.hide();
			$input.show();
		});
		this.$button.button('option', this.save_opts)
			.unbind('click')
			.click( function() {self._save();});
		this.$button_cancel.button('enable');
		this._trigger('edit');
	},
	_cancel: function()
	{
		var self = this;
		this.$button.button('option', this.edit_opts)
			.unbind('click')
			.click( function() {self._edit();});
		this.$button_cancel.button('disable');
		
		this.$form.find('.magic-item').each( function() {
			var $item = $(this);
			var $input = $item.find('.magic-input');
			var $text = $item.find('.magic-text');
			$input.hide();
			$text.show();
		});
		this._trigger('cancel');
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
	_cancel: function()
	{
		var self = this;
		//hide the errors
		this.$form.find('.magic-error').each(function() {$(this).slideUp('slow');});
		this.$form.find('.magic-item').each( function() {
			var $item = $(this);
			var $input = $item.find('.magic-input');
			var $text = $item.find('.magic-text');
			//show the text again
			$input.fadeOut('fast', function() {$text.fadeIn('fast');});
		});
		this.$button.button('option', this.edit_opts)
			.unbind('click')
			.click( function() {self._edit();});
		this.$cancel_btn.button('disable');
	},
	_data: function(data)
	{
		var self = this;
		this.$form.find('.magic-error').each(function() {$(this).slideUp('slow');});
		if(data[0] != 0) //if error
		{
			$inputs = this.$form.find('.magic-input');
			for(i in data[1].errors)
			{
				//find the input(s) to which the error relates
				$inputs.each( function() {
					var $t = $(this);
					if($t.attr('name') == i)
					{
						//show the error
						var $item = $t.closest('.magic-item').find('.magic-error');
						$item.text('' + data[1].errors[i]);
						$item.slideDown('slow');
					}
				});
			}
		}
		else
		{	
			this.$button.button('option', this.edit_opts)
				.unbind('click')
				.click( function() {self._edit();});
			this.$button_cancel.button('disable');
			this.$form.find('.magic-item').each( function() {
				var $item = $(this);
				var $input = $item.find('.magic-input');
				var $text = $item.find('.magic-text');
				//set the text
				var val = data[1].fields[$input.attr('name')];
				if(val == undefined) val = '';
				if(httppat.test(val))
				{
					val = "<a href='" + val + "'>" + val + "</a>";
					$text.html(val);
				}
				else
					$text.text(val)
				
				//show the text
				$input.hide();
				$text.show();
			});
			this._trigger('save', undefined, data);
		}
	},
});

})(jQuery);
