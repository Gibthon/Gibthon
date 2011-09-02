/*
 * jQuery UI Prompt
 *
 * Displays an alert using the jQuery UI theme
 * 
 * Written by Bill Collins for Gibthon
 */

(function ($, undefined) {

$.widget( "ui.prompt" , {
	options: {
		title: "Alert",
		message: "Body",
		modal: true,
		type: 'alert',
		confirm: {
			title: "OK",
			click: function (promptText) { return true },
			icon: "ui-icon-check"
		},
		cancel: {
			title: "Cancel",
			click: function () { return false },
			icon: "ui-icon-close"
		},
		promptText: ''
	},
	_create: function() {
		var buttons;
		var message = this.options.message;
		var self = this;
		var confirmButton = { "text":this.options.confirm.title, "click": function () { self._confirm() } };
		var cancelButton = { "text":this.options.cancel.title, "click": function () { self._cancel() } };
		console.log([confirmButton, cancelButton]);
		switch(this.options.type) {
			case 'alert':
				buttons = [ confirmButton ];
				break;
			case 'confirm':
				buttons = [ cancelButton, confirmButton ];
				break;
			case 'prompt':
				buttons = [ cancelButton, confirmButton ];
				message = message + '<br /><br /><input type="text" style="width:90%;" value="'+this.options.promptText+'" id="promptText" />';
				break;
			case 'wait':
				break;
			case 'progress':
				message = message + '<br /><br /><div id="progress" style="height:22px;"></div>'
				break;
			default:
		}
		this.element
			.html(message)
			.dialog({
				title: this.options.title,
				resizable: false,
				closeOnEscape: false,
				draggable: false,
				modal: this.options.modal,
				buttons: buttons,
				autoOpen:false,
			});
		switch($('~ div button', this.element).size()) {
			case 0:
				break;
			case 1:
				if (this.options.confirm.icon != "") {
					$('~ div button:first', this.element).button("option", "icons", {primary:this.options.confirm.icon});
				}
				break;
			case 2:
				if (this.options.confirm.icon != "") {
					$('~ div button:nth-child(2)', this.element).button("option", "icons", {primary:this.options.confirm.icon});
				}
				if (this.options.cancel.icon != "") {
					$('~ div button:first-child', this.element).button("option", "icons", {primary:this.options.cancel.icon});
				}
				break;
			default:
		}
		this.element.dialog('open');
	},
	
	_text: function () {
		return $('input', this.element).val();
	},
	
	_confirm: function () {
		if (!(this.options.type == 'prompt' && this._text() == '')) {
			this.options.confirm.click(this._text());
			this.close();
		} else {
			this._cancel();
		}
	},
	
	_cancel: function () {
		this.options.cancel.click();
		this.close();
	},

	
	close: function () {
		this.element.dialog("close");
		this._destroy();
	},
	
	_destroy: function() {
		this.element.html('');
		this.element.dialog('destroy');
		$.Widget.prototype.destroy.apply(this, arguments);
	}
});
})(jQuery);