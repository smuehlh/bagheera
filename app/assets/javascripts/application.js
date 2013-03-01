// This is a manifest file that'll be compiled into application.js, which will include all the files
// listed below.
//
// Any JavaScript/Coffee file within this directory, lib/assets/javascripts, vendor/assets/javascripts,
// or vendor/assets/javascripts of plugins, if any, can be referenced here using a relative path.
//
// It's not advisable to add code directly here, but if you do, it'll appear at the bottom of the
// the compiled file.
//
// WARNING: THE FIRST BLANK LINE MARKS THE END OF WHAT'S TO BE PROCESSED, ANY BLANK LINE SHOULD
// GO AFTER THE REQUIRES BELOW.
//
//= require jquery
//= require jquery_ujs
//= require jquery-fileupload/basic
//= require fancybox
//= require_tree .

$(function () {
    $.ajaxSetup ({
       cache: false,
       // timeout: 600000 
    });
    load_fileupload();
    load_fancy_boy();
    // $("#ajax_upload_field").validate();
});

//= require jquery.validate
//= require jquery.validate.additional-methods
// $.validator.addMethod('filesize', function(value, element, param) {
//     // param = size (en bytes) 
//     // element = element to validate (<input>)
//     // value = value of the element (file name)
//     return this.optional(element) || (element.files[0].size <= param) 
// });

function load_fileupload() {
    // $('.ajax_upload_field').validate({
    //     rules: {input:  { required: true, filesize: 1048576  }},
    //     messages: { ajax_upload_field: show_error('File must be less than 50 MB. <br /> Please contact us to upload larger files.') }
    // });
    $('.ajax_upload_field').fileupload({
        autoUpload: false,
        add: function(e, data) {
            if (data.files[0].size <= 52428800) {
                data.submit();
            }
            else {
                show_error('File must be less than 50 MB. <br /> Please contact us to upload larger files.');
            }
        },
        submit: function() {
            // first hide the predict button and results in case a subsequent submit contains false data
            $('#predict_button').hide();
            $('#options_section').hide();
            $('#results_section').hide();
            // show loading anim
            hide_show_waiting('show');
        }
    });
};

function load_fileupload_old() {
    $('.ajax_upload_field').fileupload({
        dataType: 'json',
        done: function (e, data) {
            if (data.result.error === "true") {
                show_error(data.result.message);
            }
            else {
               $('#predict_button').show();
               $('#options_section').hide();
               $('#results_section').show();
            }
        },
        submit: function() {
            // first hide the predict button in case a subsequent submit contains false data
            $('#predict_button').hide();
            $('#options_section').hide();
            $('#results_section').hide();
            // show loading anim
            hide_show_waiting('show');
        },
        always: function() {
            // hide waiting animation
            hide_show_waiting('hide');
        }
    });
};

function show_error(error) {
    $('#errors_container').html(error);
    $('#error').click();
    hide_show_waiting('hide');
};

function load_fancy_boy() {
    $("a.fancy_iframe").fancybox({
        'width'             : 530,
        'height'            : 580,
        'centerOnScroll'    : true,
        'autoScale'         : false,
        'transitionIn'      : 'elastic', //fade
        'transitionOut'     : 'elastic', //fade
        'type'              : 'iframe',
        'titlePosition'     : 'outside', //'outside', 'inside', 'over'
        'titleShow'         : false
    });

    $("a.fancy_iframe_help").fancybox({
        'width'             : 850,
        'height'            : 550,
        'centerOnScroll'    : true,
        'autoScale'         : false,
        'transitionIn'      : 'none', //fade
        'transitionOut'     : 'none', //fade
        'type'              : 'iframe',
        'titlePosition'     : 'outside', //'outside', 'inside', 'over'
        'titleShow'         : false
    });

    $("a.fancy_iframe_error").fancybox({
        'width'             : 850,
        'height'            : 550,
        'centerOnScroll'    : true,
        'autoScale'         : false,
        'transitionIn'      : 'none', //fade
        'transitionOut'     : 'none', //fade
        'type'              : 'inline',
        'titlePosition'     : 'outside', //'outside', 'inside', 'over'
        'titleShow'         : false
    });
}

function update_file_upload(ajax_upload_field, result) {
        data =  $(ajax_upload_field).data().formData;

        // result.file_id, file_name, upload_html, (protein_sequence)
        var hidden_field = data.hidden_field; //where the ID is saved
        var data_field = data.data_field; // where the upload html is displayed

        $('#'+hidden_field)[0].value = result.file_id;
        $('#'+data_field).first().html(result.upload_html);
        // show_hide_info(result.show_hide_info);
};

function hide_show_waiting(kind) {
    if (kind === 'show') {
        $('#waiting').css({'height' : $(document).height()});
        $('#waiting').show();
    }
    else {
        $('#waiting').hide();
    }
};
