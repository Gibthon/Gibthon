$(window).ready( function() {
	set_desc();
	//tabs
	$( "#tabs" ).tabs({
		ajaxOptions: {
			error: function( xhr, status, index, anchor ) {
				$( anchor.hash ).html(
					"<p>ERROR: Failed to load this tab, status: " + status + " </p>");
			}
		}
	});
	
	$('#info_btn').button({
		label: 'Edit',
		icons: {'primary': 'ui-icon-pencil'},
	}).click( function() {edit_meta();});
});

var edit_meta = function() {
	var $name = $('#c_name');
	$name.html('<input type="text" name="name" maxlength=80 size=25 value="' + $name.text() + '"/>');
	var $desc = $('#c_desc');
	$desc.html('<textarea style="margin-left:1em;" rows=5 cols=80 name="desc" >' + desc + '</textarea>');
	
	$('#info_btn').button('option', {
		'label': 'Save',
		'icons': {'primary': 'ui-icon-disk',},
	})
	.unbind('click')
	.click( function() {
		desc = $('#c_desc textarea').val();
		var name = $('#c_name input').val();
		$.getJSON('saveMeta/', {'desc': desc, 'name':name,}, function(data) {save_meta(data);});
	});	
}

var save_meta = function(data) {
	desc = $('#c_desc textarea').val();
	var name = $('#c_name input').val();
	
	$('#c_name').text(name);
	set_desc();
	
	if(data[0] != 0)
		$('#meta_error').html('<p>' + data[1] + '</p>').slideDown('slow');
	else
	{
		$('#meta_error').slideUp('slow');
		$('#status').text(data[1].modified);
	}
	
	$('#info_btn').button('option', {
		'label': 'Edit',
		'icons': {'primary': 'ui-icon-pencil',},
	})
	.unbind('click')
	.click( function() {
		edit_meta();
	});	
};

var set_desc = function() {$('#c_desc').html('<p>' + desc.replace(/\n/g, '</p><p>') + '</p>');};

