
/*
** ui.importer
* 
* Display a UI to import fragments, calling the appropriate JSON to do so
* 
*/

var importer_html = '<div id="dlg_content"></div>';

var importer_form_html = '' + 
'<div class="center-pad">' +
'	<div><button id="back_btn" /> </div>' +
'	<div id="form_holder"></div>' +
'	<div><p id="error" style="display:hidden;color:red;font-weight:bold;" class="errorlist"></p></div>' +
'	<div class="buttons" style="margin-left:auto;margin-right:auto;width:25em;text-align:center;">' +
'		<button id="cancel_btn"></button>' +
'		<button id="ok_btn"></button>' +
'	</div>' +
'</div>';

var busy = '<div style="text-align:center;margin-top:1.5em;">' +
'	<h2>Please Wait...</h2>' +
'	<img src="/static/images/spinner.gif" alt="Loading" />' +
'	<h3 id="title"></h3>' +
'</div>';

var entrez_results = '' + 
'<div class="center-pad">' +
'	<div><button id="back_btn" />  <h3 style="display:inline-block;margin-left:2.5em;">Results:</h3></div>' +
'	<div id="form_holder">' +
'	<table>' +
'		<thead>' + 
'			<tr><th>Title</th><th>Length</th>' +
' 			<th><input class="all-selected" type="checkbox" name="select-all" value="select-all"/></th></tr>' +
'		</thead>' + 
'		<tbody id="results">' +
'		</tbody>' +  
'	</table></div>' +
'	<div><p id="error" style="display:hidden;color:red;font-weight:bold;" class="errorlist"></p></div>' +
'	<div class="buttons" style="margin-left:auto;margin-right:auto;width:25em;text-align:center;">' +
'		<button id="cancel_btn"></button>' +
'		<button id="ok_btn"></button>' +
'	</div>' +
'</div>';

(function( $, undefined ) {

$.widget("ui.importer", {
	_create: function() {
		var self = this;
		this.updated = false; //set to true to reload the page on close		
		
		this.$dlg = $(this.element[0])
			.html(importer_html)
			.dialog({
				autoOpen:false,
				resizable:false,
				modal:true,
				title:"Add new Fragment",
				close:function () {
					self._show_main();
					if(self.updated)
						location.reload();
				},
				height:"260",
				width:"640"
			});
		this.$content = this.$dlg.find('#dlg_content');
		this._show_main();
		this._results = new Array();
	},
	open: function(){ this.$dlg.dialog('open');},
	close: function() { this.$dlg.dialog('close');},		
	_show_main: function() { //show the main window
		var self = this;
		var $new = $(document.createElement('div'));
		$new.load('/fragment/import/', {}, function () {
			self._normal_size();
			$new.find('#parts').click(function() {
				self._show_parts();
			});
			$new.find('#entrez').click(function() {
				self._show_entrez();
			});
			$new.find('#upload').click(function() {
				self._show_upload();
			});
			$new.find('#manual').click(function() {
				self._show_manual();
			});
			self._show($new);
		});
	},
	_show_busy: function(title)
	{
		if(title == undefined)
			title = "";
		var $busy = $(document.createElement('div'));
		$busy.html(busy);
		$busy.find('#title').text(title);
		this.$content.html($busy);
	},
	_show_parts: function(){
		var self = this;
		var $new = $(document.createElement('div')).html(importer_form_html);
		$new.find('#form_holder').load('/fragment/import/part/', {}, function(){
			$new.find('#add_form').submit(function() {
				return false;
			});
			$new.find('#extender').formExtender();
			self._show($new);
			self._auto_size();
		});
		
		$new.find('#ok_btn').button({
				label: 'Import',
			icons: {primary: 'ui-icon-gear'},
		}).click(function() {
			self.updated = true;
			$new.find('#extender-copy input').each(function() {
				var name = $(this).val();
				var $state = $(this).siblings('.status');
				if($state.hasClass('status-ok')) //we've already imported this one!
					return; //let's not do it again!
				var $error = $(this).siblings('.extender-error').slideUp('fast');
				
				$state
					.removeClass('status-error')
					.addClass('status-busy');
				$.getJSON('import/part/go/', {'part': name,}, function(data){
					console.log("DATA: status:" + data[0] + " value:" + data[1])
					if(data[0] != 0) //if error
					{
						$state.removeClass('status-busy').addClass('status-error');
						console.log("Error: " + data[1]);
						$error.text(data[1]).slideDown('slow');
					}
					else
					{
						$state.removeClass('status-busy').addClass('status-ok');
					}
				}); 
				
			});
		});
		this._enable_back_cancel($new, 'Close');
		
	},
	_show_manual: function(){
		var self = this;
		var $new = $(document.createElement('div')).html(importer_form_html);
		$new.find('#form_holder').load('/fragment/import/manual/', {}, function(){
			var $form = $new.find('#add_form').submit(function() {
				//submit via AJAX instead
				$.getJSON($form.attr('action'), $form.serialize(), function(data) {
					if(data[0] != 0) //error
					{
						console.log("Error: " + data[1]);
						for( key in data[1])
						{
							error = data[1][key];
							console.log(key + ':' + error);
						console.log("$new.find('#" + key + "_error').length: " + $new.find('#' + key +'_error').length);
							$new.find('#' + key +'_error').text('ERROR: ' + error).slideDown('fast');
						}
					}
					else//success
					{
						location.reload();
					}
				});
				return false;
			});
			self._show($new);
			self._auto_size('top');
		});
		
		$new.find('#ok_btn').button({
				label: 'Save',
			icons: {primary: 'ui-icon-disk'},
		}).click( function() {
			var error = false;
			var $name = $new.find('#id_name');
			if(($name.val() == ''))
			{
				error = true;
				$new.find('#name_error').text('Name cannot be empty.').slideDown('fast');
			}
			else
				$new.find('#name_error').slideUp('fast');
			var $desc = $new.find('#id_desc');
			var $seq = $new.find('#id_seq');
			if(($seq.val() == ''))
			{
				error = true;
				$new.find('#seq_error').text('Sequence cannot be empty.').slideDown('fast');
			}
			else
				$new.find('#seq_error').slideUp('fast');
			
			if(error) return;
			
			console.log('submit');
			$new.find('form').submit();
		});
		
		this._enable_back_cancel($new);
	},
	_show_entrez: function(error){
		var self = this;
		var $new = $(document.createElement('div')).html(importer_form_html);
		
		$new.find('#form_holder').load('/fragment/import/entrez/', {}, function(){
			$new.find('#add_form').submit(function() {
				self._entrez_search();
				return false;
			});
			//show the whole kaboodle
			self._show($new);
		});
		if((error != undefined) && (error != ""))
		{
			$new.find('#error')
				.text(error)
				.css('display', 'block');
		}
		//hookup the buttons
		$new.find('#ok_btn').button({
			label: 'Search',
			icons: {primary: 'ui-icon-search'},
		}).click(function() { 
			self._entrez_search();
		});
		this._enable_back_cancel($new);
		
	},
	_show_upload: function(){
		
		var self = this;
		var $new = $(document.createElement('div')).html(importer_form_html);
		
		$new.find('#form_holder').load('/fragment/import/upload/', {}, function(){
			$new.find('#add_form').submit(function() {
				return false;
			});
			//show the whole kaboodle
			self.$content.fadeOut('fast', function() {
				$(this).html($new)
					.find('#fileupload').fileupload();
				$(this).fadeIn('fast');
				self._auto_size();
			});
			self.updated = true;
		});
		
		//hookup the buttons
		$new.find('#ok_btn').button({
			label: 'Done',
			icons: {primary: 'ui-icon-close'},
		}).click(function() { 
			location.reload();
		});
		this._enable_back_cancel($new);
		
	},
	
	_enable_back_cancel: function($new, cancel_label){
		if (cancel_label == undefined) cancel_label = 'Cancel';
		var self = this;
		$new.find('#back_btn').button({
			label: 'Back',
			icons: {primary: 'ui-icon-arrowreturn-1-w'},
		}).click(function() { self._show_main();});
		$new.find('#cancel_btn').button({
			label: cancel_label,
			icons: {primary: 'ui-icon-cancel'},
		}).click(function() { self.close(); });
	},
	_show: function($new){
		this.$content.fadeOut('fast', function() {
			$(this)
				.html($new)
				.fadeIn('fast');
		});
	},

//######################################################## ENTREZ
	_entrez_search: function(){
		var self = this;
		var data = this.$content.find('#add_form').serialize();
		this._show_busy("Searching NCBI Entrez...");
		$.getJSON('/fragment/import/entrez/search/', data, function(data) {
			if(data[0] != 0)
				self._show_entrez(data[1]);
			else
				self._entrez_results(data[1]);
		});
	},
	_entrez_results: function(ids){
		var self = this;
		var cmds = new Array();
		for( i in ids)
		{
			cmds.push({'desc': "Fetching info for id '" + ids[i] + "'", 
				'url': '/fragment/import/entrez/summary/', 
				'data': {'id': ids[i],},});
		} 
		var s = '';
		if(ids.length > 1) s = 's';
		this.$content.loader({
			'commands': cmds, 
			'autoStart': true,
			'data': function(val) {self._entrez_summary(val);},
			'done': function(val) {self._entrez_show_results(val);},
			'title': 'Found ' + ids.length + ' result' + s + ', fetching summary information...',
		});
	},
	_entrez_summary: function(summary){
		this._results.push(summary);
	},
	_entrez_show_results: function(errors){//ignore errors and show results
		var self = this;
		var $results = $(document.createElement('div'));
		$results.html(entrez_results);
		var $results_tbl = $results.find('#results');
		
		for (i in this._results)
		{
			var s = this._results[i];
			$results_tbl.append(make_row(s));
		}
		
		this._enable_back_cancel($results);
		var $ok = $results.find('#ok_btn').button({
			label: 'Import',
			icons: {primary: 'ui-icon-transfer-e-w'}, //ui-icon-gear
			disabled: true,
		});
			
		$fl = $results.find('table')
			.fragmentList({
				selectChanged: function(event, d){
					if(d.selected == 0)
						$ok.button('disable');
					else
						$ok.button('enable');
				}	
			})
			.dataTable({
				"bAutoWidth": false,
				"bJQueryUI": true,
				"aoColumns": [
					{
						"sWidth":"674px"
					},
					{
						"sWidth":"100px",
						"sType": 'numeric'
					},
					{
						"sWidth":"40px"
					},
				]
			});		
		$ok.click(function() { 
			//if some fragments have been selected, import them
			var ids = $fl.fragmentList('getList');
			var commands = new Array()
			for( i in ids )
			{
				var id = ids[i];
				commands.push( {'desc': "Importing id '" + id + "'...", 
								'url':'/fragment/import/entrez/import/',
								'data': {'id': id,}, });
			}
			var s = '';
			if(ids.length > 1) s = 's';
			self.$content.loader({
				'commands': commands, 
				'autoStart': true,
				'data': function(val) {},
				'done': function(val) {location.reload();},
				'title': 'Importing ' + ids.length + ' fragment' + s + '...',
			});
			self._normal_size();
		});	
		
		this._auto_size();
		this.$content.html($results);
	},
	
	
//DIALOG RESIZING
	_full_size: function()
	{
		//set position and height;
		this.$dlg.dialog('option', 'height',  ($(window).height() - 10));
		this.$dlg.dialog('option', 'position','center');
	},
	_normal_size: function()
	{
		//set position and height;
		this.$dlg.dialog('option', 'height',  260);
		this.$dlg.dialog('option', 'position','center');
	},
	_auto_size: function(position)
	{
		if(position == undefined) position = 'center';
		this.$dlg.dialog('option', 'height',  'auto');
		this.$dlg.dialog('option', 'position',position);
	},
});

})(jQuery);

var make_row = function(summary)
{
	return '' +
	'<tr>' +
	'	<td>' +
	'		<p class="table-para" >' + summary.Title + '</p>' +
	'	</td>' +
	'	<td>' +
	'		<p class="table-para" >' + summary.Length + '</p>' +
	'	</td>' +
	'	<td>' +
	'		<p class="table-para">' +
	'		<input class="selected-check" id="' + summary.Id + '-checkbox" ' +
	'				type="checkbox" name="selected" value="' + summary.Id + '" />' +
	'		</p>' +
	'	</td>' +
	'</tr>';
}

/*
 * ui.loader
 * 
 * 	Performs a list of AJAX actions while displaying a progress bar
 *  Displays the error message if one accurs.
 *  Actions which return values (other than OK or ERROR) are not supported
 * 
 *  usage: $(SELECTOR).loader({commands: [{	'desc': 'Command Description',
 * 											'url': '/url/to/call/',
 * 											'data': {	'the': 'data to post',
 * 														'to': 'the url'}, 
 * 											}, MORE COMMANDS HERE];
 * 
 * then call .loader('start') to begin, or set option autoStart to true.
 */

var loader_html = '' + 
'<div style="text-align:center;margin-top:1.5em;">' +
'	<h2>Please Wait...</h2>' +
'	<img src="/static/images/spinner.gif" alt="Loading" style="margin-bottom:.5em;" />' +
'	<h3 id="title"></h3>' +
'	<h4 id="current"></h4>' +
'	<h5 id="text_progress" style="margin-bottom:.2em;"></h5>' +
'	<div id="errors" style="color:red;"></div>' +
'	<div id="loader_progress" style="height:1.1em;"></div>' +
'	<div id="buttons" style="display:none;"><button id="close_btn"></button></div>' +
'</div>';

(function( $, undefined ) {

$.widget("ui.loader", {
	options: {
		commands: [], //a list of commands, in the form [{'desc': DESCRIPTION, 'url': URL, 'data': DATA}, etc]
		autostart: false, //should the actions start immediately?
		data: function(data) {}, //what to do with the data
		done: function(errors) {}, //errors is the number of errors
		title: '',
	},
	_init: function()
	{
		if(this.options.commands == [])
		{
			console.error('ui.loader: init called without commands!');
			return;
		}
		
		var self = this;
		var $new = $(document.createElement('div')).html(loader_html);
		this.$el = $(this.element[0]).fadeOut('fast', function() {
			$(this).html($new).fadeIn('fast', function() {
				if(self.options.autoStart)
				{
					self._started = true;
					self._next_action();
				}
			});
		});
		this.$progress = $new.find('#loader_progress').progressbar({value: 0,});
		this.$text_progress = $new.find('#text_progress').text('0 / 0');
		this._cmd = this.options.commands;
		this._len = this.options.commands.length;
		this._next = 0; //the offset of the function to be performed next
		this._errors = 0;
		this.$title = $new.find('#title').html(this.options.title);
		this.$current = $new.find('#current');
		this.$errors = $new.find('#errros');
		
		this._started = false;
	},
	start: function() //start performing the commands
	{
		if(!this._started)
		{
			this._started = true;
			this._next_action();
		}
	},
	progress: function() //return the progress
	{
		return this._next;
	},
	_next_action: function()
	{
		//perform the next action
		if(this._next >= this._len)
		{
			this._on_complete();
			return;
		}
		var self = this;
		var cmd = this._cmd[this._next];
		this.$current.text(cmd.desc);
		var type = cmd.type;
		if(type == undefined)
			type = 'GET'; //default to GET
		$.ajax({
			'url': cmd.url,
			'data': cmd.data,
			'type': type,
			'success': 
				function(data){
					if(data[0] != 0)
					{
						self.$errors.append(data[1] + '<br/>');
						this._errors = this._errors + 1;
					}
					
					//move on to the next command
					self._next = self._next + 1;
					self._update_progressbar();
					self._next_action();
					
					if(data[0] == 0)
						self.options.data(data[1]); //call data callback
				},
		});
	},
	_on_complete: function()
	{
		this.options.done(this._errors);		
	},
	_update_progressbar: function()
	{
		this.$text_progress.text(this._next + ' / ' + this._len);
		this.$progress.progressbar('value', 100 * (this._next / this._len));
	},
});

})(jQuery);
