function makeHttpObject() {
  try {return new XMLHttpRequest();}
  catch (error) {}
  try {return new ActiveXObject("Msxml2.XMLHTTP");}
  catch (error) {}
  try {return new ActiveXObject("Microsoft.XMLHTTP");}
  catch (error) {}

  throw new Error("Could not create HTTP request object.");
}

var validate = function(input) {
	re = /^\d+(?:\.\d*)?$/;
	if(re.exec(input.value) === null){
		$(input.parentNode).addClass("errorcell");
		return false;
	}else{
		$(input.parentNode).removeClass("errorcell");
		return true;
	}
}

var Mix = function () {
	this.enz_d = $('#enzyme_desired')[0];
	this.enz_s = $('#enzyme_stock')[0];
	this.buff_d = $('#buffer_desired')[0];
	this.buff_s = $('#buffer_stock')[0];
	this.dntp_d = $('#dntp_desired')[0];
	this.dntp_s = $('#dntp_stock')[0];
	this.prim_d = $('#primer_desired')[0];
	this.temp_d = $('#template_desired')[0];
	this.repeats = $('#repeats')[0];
	this.error = $('#error_margin')[0];
	this.vol_e = $('#volume_each')[0];
	var mixform = $('form.mixform');
	this.errors = 0;
	this.warnings = 0;
	this.fragments = []
	for (f in mixform){
		if (mixform[f].constructor == HTMLFormElement){
			this.fragments.push(new MixFragment(this, mixform[f]));
		}
	}
	this.go();
}

Mix.prototype.validate = function () {
	this.errors = 0;
	for (x in this){
		if (this[x].constructor == HTMLInputElement){
			this.errors += !validate(this[x]);
		}
	}
}

Mix.prototype.m = function () {
	return this.repeats.value * (1 + this.error.value/100);
}

Mix.prototype.go = function () {
	this.warnings = 0;
	this.validate();
	for (f in this.fragments) {
		this.fragments[f].validate();
	}
	if (this.errors > 0) {
		return false;
	}
	for (f in this.fragments) {
		this.fragments[f].go();
		this.fragments[f].warn();
	}
	if (this.warnings > 0) {
		$('#warning').show();
	} else {
		$('#warning').hide();
	}
}


var MixFragment = function (mix, form) {
	this.mix = mix;
	this.prim_t_s = form.elements['primer_top_stock'];
	this.prim_b_s = form.elements['primer_bottom_stock'];
	this.temp_s = form.elements['template_stock'];
	this.enz_v = form.elements['enzyme_vol'];
	this.buff_v = form.elements['buffer_vol'];
	this.dntp_v = form.elements['dntp_vol'];
	this.prim_t_v = form.elements['primer_top_vol'];
	this.prim_b_v = form.elements['primer_bottom_vol'];
	this.temp_v = form.elements['template_vol'];
	this.water_v = form.elements['water_vol'];
	this.total_v = form.elements['total_vol'];
}

MixFragment.prototype.go = function () {
	this.enz_v.value = (this.mix.m()*this.mix.enz_d.value/this.mix.enz_s.value).toFixed(1);
	this.buff_v.value = (this.mix.vol_e.value*this.mix.m()*this.mix.buff_d.value/this.mix.buff_s.value).toFixed(1);
	this.dntp_v.value = (this.mix.vol_e.value*this.mix.m()*this.mix.dntp_d.value/this.mix.dntp_s.value).toFixed(1);
	this.temp_v.value = (this.mix.m()*this.mix.temp_d.value/this.temp_s.value).toFixed(1);
	this.prim_t_v.value = (this.mix.vol_e.value*this.mix.m()*this.mix.prim_d.value/this.prim_t_s.value).toFixed(1);
	this.prim_b_v.value = (this.mix.vol_e.value*this.mix.m()*this.mix.prim_d.value/this.prim_b_s.value).toFixed(1);
	this.water_v.value = (this.mix.m()*this.mix.vol_e.value - (parseFloat(this.enz_v.value) + parseFloat(this.buff_v.value) + parseFloat(this.dntp_v.value) + parseFloat(this.temp_v.value) + parseFloat(this.prim_t_v.value) + parseFloat(this.prim_b_v.value))).toFixed(1);
	this.total_v.value = (this.mix.m()*this.mix.vol_e.value).toFixed(1);
}

MixFragment.prototype.validate = function () {
	for (x in this){
		if (this[x].constructor == HTMLInputElement && !this[x].readOnly){
			this.mix.errors += !validate(this[x]);
		}
	}
}

MixFragment.prototype.warn = function () {
	for (x in this){
		if (this[x].constructor == HTMLInputElement && this[x].readOnly){
			if (parseFloat(this[x].value) < 1 ){
				$(this[x].parentNode).addClass("warningcell");
				this.mix.warnings += 1;
			} else {
				$(this[x].parentNode).removeClass("warningcell");
			}
		}
	}
}

$('document').ready(function () {
	$('#protocol-accordion').accordion({
		collapsible:true,
	});
	$('#warning').hide();
	$('.mix').change(function(){
		validate(this);
		if ($('.errorcell').size() > 0){
			$('button#download_protocol').button('disable');
		} else {
			$('button#download_protocol').button('enable');
		}
	});
	$('button#reset_offset').button({
		icons:{primary:'ui-icon-refresh'}
	}).click(function () {
		$('#primer_offset_confirm').html('You are about to reset all of the primers to 0 offset.');
		$('#primer_offset_confirm').dialog( "option", "buttons", { 
			"Ok": function() {
				$(this).dialog("close"); 
				$('#process').dialog('open');
				var h = makeHttpObject();
				h.open('GET', 'reset', true)
				h.onreadystatechange = function (){
					p = parseInt(this.responseText.split(':').pop());
					$('#process_progress').progressbar('value', p);
					if (p==100){
						window.location.reload();
					}
				}
				h.send(null);
			},
			"Cancel": function() { 
				$(this).dialog("close"); 
			} 
		});
		$('#primer_offset_confirm').dialog('open');
	});
	$('.primer_info_wrapper').hide();
	$('#process_progress').progressbar({ value: 0 });
	$('#process').dialog({
		autoOpen:false,
		modal:true,
		resizable:false,
		title:'Please wait',
		closeOnEscape:false,
		draggable:false,
		open: function(event, ui) {
			$(".ui-dialog-titlebar-close").hide();
		},
		close: function(event, ui) {
			$('#process_progress').progressbar('value', 0);
		}
	});
	$('#primer_offset_confirm').dialog({
		autoOpen:false,
		width:"500",
		modal:true,
		resizable:false,
		title:'Confirm offset',
		closeOnEscape:false,
		draggable:false,
	});
	$('#wait').dialog({
		autoOpen:false,
		modal:true,
		resizable:false,
		title:'Please wait',
		closeOnEscape:false,
		draggable:false,
		open: function(event, ui) {
			$(".ui-dialog-titlebar-close").hide();
		},
		close: function(event, ui) {
			$('.ui-dialog-titlebar-close').show();
		}
	});
	
	$('button.primer_top').button({
		icons:{primary:'ui-icon-circle-arrow-w'}
	});
	$('button.primer_bottom').button({
		icons:{secondary:'ui-icon-circle-arrow-e'}
	});
});