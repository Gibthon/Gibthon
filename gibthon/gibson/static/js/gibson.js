/* prevents collisions between accordion and sortable animations */
var stop = false;

/* save function */
var save = function () {
	var feature_select = [];
	var direction = [];
	$('.formthing').each(function(i,e){
		feature_select.push([ e.elements[0].value, e.elements[1].value ]);
		direction.push(e.elements[2].checked ? e.elements[2].value : e.elements[3].value);
	});
	var order = [];
	$('#fragment_list > div').each(function(i,e){order.push(e.id)});
	$('#status').load('save',{
		'order[]':order,
		'feature_select[]':feature_select,
		'direction[]':direction
	});
	$('#summary').load('summary');
};

/* http factory for process request */
function makeHttpObject() {
  try {return new XMLHttpRequest();}
  catch (error) {}
  try {return new ActiveXObject("Msxml2.XMLHTTP");}
  catch (error) {}
  try {return new ActiveXObject("Microsoft.XMLHTTP");}
  catch (error) {}

  throw new Error("Could not create HTTP request object.");
}

var refresh = function () {
	$('#fragment_list').hide().accordion('destroy');
	$('#fragment_list').load('fragments', function() {
		$(this).show();
		$('button[id^="delete/"]')
			.button({ icons:{primary:'ui-icon-trash'} })
			.click(function(event){
				var targetUrl = this.id;
				$('#fragment_delete').dialog({
					buttons : {
						"Cancel" : function() {
							$(this).dialog("close");
						},
						"Confirm" : function() {
							$.get(targetUrl, function () { refresh(); $('#fragment_delete').dialog('close'); });
						}
					}
				});
				$('#fragment_delete').dialog('open');
			});
		$('.formthing').each(function(i,form) {
			form.elements[2].parentNode.id = form.elements[2].name;
			$(form.elements[2].parentNode).buttonset();
			if (form.elements[0].value != form.elements[1].value)
			{
				$(form.elements[2].parentNode).buttonset('disable');
			}
		});
		$('input[type="radio"]').change(function() {
			save();
		});
		$('select').change(function() {
			var e = this.form.elements;
			if (e[0].selectedIndex != e[1].selectedIndex)
			{
				$(e[2].parentNode).buttonset('disable');
				if (e[0].selectedIndex < e[1].selectedIndex)
				{
					e[2].checked = true;
					e[3].checked = false;
					$(e[2].parentNode).buttonset();
				} else {
					e[2].checked = false;
					e[3].checked = true;
					$(e[2].parentNode).buttonset();
				}
			} else {
				$(e[2].parentNode).buttonset('enable');
			}
			save();
		});
		$('button[id^="view"]')
			.button({ icons:{primary:'ui-icon-search'} })
			.click(function(){
				$('#fragment_viewer_content').load(this.id);
				$('#fragment_viewer').dialog('open');
			});
			
		/* fragment accordion */
		$('#fragment_list > div > h3').click(function(event){
			if(stop) {
				event.stopImmediatePropagation();
				event.preventDefault();
				stop = false;
			}
		});
		$('#fragment_list')
			.accordion({
				collapsible:true,
				header: "> div > h3"
			})
			.sortable({
				axis: "y",
				handle: "h3",
				stop: function() {
					stop = true;
					save()
				}
			});
		$('#summary').load('summary');
	});
}
	

$(document).ready(function() {
	/* load the fragments */
	refresh();
	
	/* general initialisations */
	$('#process_progress').progressbar({ value: 0 });
	
	/* load the settings */
	$('div#settings-wrapper').load('settings');
	
	/* buttons */
	$('button#close_progress')
		.button({ icons:{primary:'ui-icon-close'} })
		.click(function(){ $('#wait').dialog('close'); });
	$('button#primers_progress')
		.button({ icons:{primary:'ui-icon-transferthick-e-w'} })
		.click(function(){ window.location.href='primers'; });
	$('button#info')
		.button({ icons:{primary:'ui-icon-pencil'} })
		.click(function(){ $('#construct_edit').dialog('open'); });
	$('button#download')
		.button({ icons:{primary:'ui-icon-link'} })
		.click(function(event){ window.location.href = '../{{ construct.name }}.gb'; });
	$('button#info_edit')
		.button({ icons:{primary:'ui-icon-check'} })
		.click(function(event){
			$.post('', $('form[action^="info_edit"]').serialize(), function(url){
				window.location.href = url;
			});
		});
	$('button#primers')
		.button({ icons:{primary:'ui-icon-transferthick-e-w'} })
		.click(function(){ window.location.href="primers" });
	$('button#process')
		.button({ icons:{primary:'ui-icon-check'} })
		.click( function(event) {
			save();
			$('#wait').dialog('open');
			var h = makeHttpObject();
			h.open('GET', 'process', true)
			h.onreadystatechange = function (){
				p = parseInt(this.responseText.split(':').pop());
				$('#process_progress').progressbar('value', p);
				if (p==100){
					$('button.button_progress').button('enable');
				}
			}
			h.send();			
		});
	$('button#delete')
		.click(function(){ $('#construct_delete').dialog('open'); })
		.button({ icons:{primary:'ui-icon-trash'} });
	$('button#save')
		.button({ icons:{primary:'ui-icon-disk'} })
		.click(function(){ save() });
	$('button#add')
		.button({ icons:{primary:'ui-icon-plusthick'} })
		.click(function(event){
			$('#fragment_browser_content').load('add');
			$('#fragment_browser').dialog('open');
		});
		
	$('button.button_progress').button('disable');
		
	
	/* dialogs */
	$('#construct_edit').dialog({
		autoOpen:false,
		resizable:false,
		modal:true,
		title:'Add new Construct',
		width:400
	});
	$('#wait').dialog({
		autoOpen:false,
		modal:true,
		resizable:false,
		title:'Please wait',
		closeOnEscape:false,
		draggable:false,
		open: function(event, ui) {
			$('button.button_progress').button('disable');
			$(".ui-dialog-titlebar-close").hide();
		},
		close: function(event, ui) {
			$('#process_progress').progressbar('value', 0);
			$('button.button_progress').button('disable');
		}
	});
	$('#fragment_delete').dialog({
		autoOpen:false,
		modal:true,
		resizable:false,
		title:'Confirm Delete'
	});
	$('#construct_delete').dialog({
		autoOpen:false,
		modal:true,
		resizable:false,
		title:'Confirm Delete',
		buttons : {
			"Cancel" : function() {
				$(this).dialog('close');
			},
			"Delete" : function() {
				window.location.href = 'delete';
			}
		}
	});
	$('#fragment_browser').dialog({
		autoOpen:false,
		width:600,
		modal:true,
		resizable:false,
		title:'Add Fragment'
	});
	
	$('#fragment_viewer').dialog({
		autoOpen:false,
		draggable:false,
		height:600,
		width:1000,
		modal:true,
		resizable:false,
		title:'Browse Fragment',
		close:function(){
			$('#fragment_viewer_content').html('');
		}
	});
	
	
});