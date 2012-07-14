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
		this.$copy = this.$el.find('#extender-copy').detach();
		this.$target = this.$el.find('.extender-target');
		
		this.$el.find('button.extender-add').button({
			label: 'Add',
			icons: {primary: 'ui-icon-plusthick',},
		}).click( function(event) {
			self._add(event);
			return false;
		});
		
		this.$el.find('button.extender-remove').each( function() {
			$(this).button({
				text: false,
				icons: {primary: 'ui-icon-trash',},
			}).click( function(event) {
				self._remove(event);
				return false;
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
	Add: function(values) {this._add(undefined, values);},
	_add: function(event, values){
		var self=this;
		var $clone = this.$copy.clone();
		$clone.addClass('extender-item').attr('id', '');
		
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
			return false;
		});
		
		if(values != undefined)
		{
			$clone.find('.magic-input').each( function(i, v) {
				if(i < values.length)
				{
					$(v).val(values[i]);
				}
			});
		}
		
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
		submit_on_save: true, //whether to submit the form on save, 
        //  or whether to allow the save handler to do it
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
		this.state = 'default';
		if((this.method == '') || (this.method == undefined))
			this.method = 'GET'
		
		this.edit_opts = {'label': 'Edit', 'icons': {'primary': 'ui-icon-pencil',},}
		this.save_opts = {'label': 'Save', 'icons': {'primary': 'ui-icon-disk',},}
		this.$button = this.$form.find('button.magic-button')
			.button(this.edit_opts)
			.click( function() {return self._edit();});
		this.$button_cancel = this.$form.find('button.magic-button-cancel')
			.button({label: 'Cancel', icons: {primary: 'ui-icon-cancel'}, disabled: true,})
			.click( function() {return self._cancel();});
			
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
		this.state = 'editable';
		this.$form.find('.magic-item').each( function() {
			self._set_editable($(this));
			
		});
		this.$button.button('option', this.save_opts)
			.unbind('click')
			.click( function() {return self._save();});
		this.$button_cancel.button('enable');
		this._trigger('edit');
		return false;
	},
	_save: function()
	{
		var self = this;
		self.state = 'default';
		if(this.options.submit_on_save)
		{
			$.ajax({
				url: this.action,
				dataType: 'json',
				data: this.$form.serialize(),
				type: this.method,
				success: function(data) {self._data(data);},
			});
		}
		else
		{
			self._data([0, undefined]);
		}
		return false;
	},
	_cancel: function()
	{
		var self = this;
		//hide the errors
		self.state = 'default';
		this.$form.find('.magic-item').each( function() {
			self._set_default($(this));
			
		});
		this.$button.button('option', this.edit_opts)
			.unbind('click')
			.click( function() {self._edit();});
		this.$button_cancel.button('disable');
		this._trigger('cancel');
		return false;
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
				.click( function() {self._edit(); return false;});
			this.$button_cancel.button('disable');
			this.$form.find('.magic-item').each( function() {
				var $item = $(this);
				var $input = $item.find('.magic-input');
				
				//set the text
				var val = '';
				if( data[1] != undefined)
				{
					val = data[1].fields[$input.attr('name')];
					if(val == undefined) val = '';
					self._set_default($item, val);
				}
				else 
					self._set_default($item);
				
			});
			this._trigger('save', undefined, data);
		}
		return false;
	},
	SetState: function($item)
	{
		if(!$item.hasClass('magic-item'))
		{
			console.log('magicForm.SetState given bogus item');
			return;
		}
		switch(this.state)
		{
			case 'default':
				this._set_default($item);
				break;
			case 'editable':
				this._set_editable($item);
				break;
			default:
				return;
		}
	},
	_set_editable: function($item)
	{
		var $input = $item.find('.magic-input');
		var $text = $item.find('.magic-text');
		//set the value
		$input.val($text.text());
		$text.hide();
		$input.show();
	},
	_set_default: function($item, value)
	{
		var $input = $item.find('.magic-input');
		var $text = $item.find('.magic-text');
		
		
		if(value != undefined)
		{
			if(httppat.test(value))
			{
				val = "<a href='" + value + "'>" + value + "</a>";
				$text.html(value);
			}
			else
				$text.text(value)
		}
		else
			$text.text($input.val());
		
		$item.find('.magic-error').slideUp('slow');
		//show the text again
		$input.hide();
		$text.show();
	},
});

})(jQuery);
