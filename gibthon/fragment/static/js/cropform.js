/*
 * Widget to control the cropping form
 *
 */

$.widget('ui.cropform', {
    options: {
        start: 0,
        end: 0,
        fragment: null,
        cancel: function() {},
    },
    _create: function()
    {
        var self = this;
        var el = $(this.element[0]);
        this.el = el;
        this.new_settings = el.find('#new_settings');
        this.result_type = el.find('select#result_type');
        this.cb_internal = el.find('input#cb_internal');
        this.cb_all = el.find('input#cb_all');
        this.disable = el.find('div#disableme');

        el.find('#cr_cancel_btn').button({
            label: 'Cancel',
            icons: {primary: 'ui-icon-cancel',},
        }).click(function() {
            self._trigger('cancel');
        });
        el.find('#cr_crop_btn').button({
            label: 'Crop',
            icons: {primary: 'ui-icon-scissors',},
        }).click(function() {self._crop();});

        this.result_type.selectmenu({
            style:'dropdown',
            width: 200,
            change: function(ev, ui)
            {
                var act = $(ui.option).val();
                if(act == 'new')
                    self.enableNewSettings();
                else if(act == 'overwrite')
                    self.disableNewSettings();
            },
        });

        this.cb_internal.change(function(){
            if(self.cb_internal.is(':checked'))
                self.disable.slideDown(200);
            else
                self.disable.slideUp(200);
        });

        var o = this.options;
        this.name = $('#new_name').val(o.fragment.getName()+'_crop'); 
        this.desc = $('#new_desc').val(o.fragment.getDesc()); 
        this.start = $('#range_from').text(this.options.start);
        this.end = $('#range_to').text(this.options.end);

    },
    setRange: function(start, end)
    {
        this.start.text(start);
        this.end.text(end);
        this._updateDesc();
    },
    _updateDesc: function()
    {
        this.desc.val(this.options.fragment.getDesc() + '. Cropped to '+
                      this.start.text()+'bp -> '+this.end.text()+'bp.');
    },
    disableNewSettings: function()
    {
        this.new_settings.slideUp(500);
    },
    enableNewSettings: function()
    {
        this.new_settings.slideDown(500);
    },
    _crop: function()
    {
        this.el.find('#beforecrop').slideUp(500);
        this.el.find('#cropwait').slideDown(500);

        this.options.fragment.crop({
            start: this.start.text()-1,
            end: this.end.text()-1,
            f_internal: this.cb_internal.is(':checked') ? 1 : 0,
            f_all: (this.cb_internal.is(':checked') && 
                    this.cb_all.is(':checked')) ? 1 : 0,
            result: this.result_type.selectmenu('value'),
            new_name: this.name.val(),
            new_desc: this.desc.val(),
        }, function(id){
            if($('#go_new').is(':checked'))
                window.location = '/fragment/'+id;
            else
                location.reload();
        });

    }
});

