function pathname = uigetdir_java(varargin)

%import java.awt.*;
%import java.awt.event.*;
%import javax.swing.*;
%import javax.swing.filechooser.*;
%import java.lang.String

if nargin == 0 
  initial_dir = java.lang.String(pwd);
  dlgtitle = 'Select directory';
elseif nargin >= 1 & exist(varargin{1},'dir')
  initial_dir = java.lang.String(varargin{1});
  dlgtitle = 'Select directory';
else
  errordlg('Input argument must be a valid directory',...
      'Input Argument Error!')
  return
end
if nargin >= 2 & ischar(varargin{2})
  dlgtitle = varargin{2};
end

initial_dirf = java.io.File(initial_dir);
TRUE = 1;
FALSE = 0;

% javax.swing.UIManager.setLookAndFeel('com.sun.java.swing.plaf.windows.WindowsLookAndFeel');
javax.swing.UIManager.setLookAndFeel('javax.swing.plaf.metal.MetalLookAndFeel');
% See note below for more setable options
javax.swing.UIManager.put('FileChooser.fileNameLabelText', 'Directory name');

f=java.awt.Frame;
j=javax.swing.JFileChooser(initial_dir);
j.setFileSelectionMode(j.DIRECTORIES_ONLY);
j.setFileHidingEnabled(TRUE); 
j.setCurrentDirectory(initial_dirf);
j.setSelectedFile(initial_dirf);

str = 'Select';
jstr = java.lang.String(str);
j.setDialogTitle(dlgtitle);
j.setApproveButtonToolTipText(dlgtitle);
j.setApproveButtonText('Select');
j.setApproveButtonMnemonic('S');

returnVal = j.showDialog(f,jstr)
if returnVal == javax.swing.JFileChooser.APPROVE_OPTION
  pathnamef = j.getSelectedFile;
  pathname = char(pathnamef.getAbsolutePath)
else
  pathname = [];
  disp('Open command cancelled by user.')
end

%      Can also set Filechooser properties for:
%                     FileChooser.lookInLabelText
%                     FileChooser.upFolderToolTipText
%                     FileChooser.filesOfTypeLabelText
%                     FileChooser.fileNameLabelText
%                     FileChooser.homeFolderToolTipText
%                     FileChooser.newFolderToolTipText
%                     FileChooser.listViewButtonToolTipTextlist
%                     FileChooser.detailsViewButtonToolTipText
%                     FileChooser.saveButtonText=Save
%                     FileChooser.openButtonText=Open
%                     FileChooser.cancelButtonText=Cancel
%                     FileChooser.updateButtonText=Update
%                     FileChooser.helpButtonText=Help
%                     FileChooser.saveButtonToolTipText=Save
%                     FileChooser.openButtonToolTipText=Open
%                     FileChooser.cancelButtonToolTipText=Cancel
%                     FileChooser.updateButtonToolTipText=Update
%                     FileChooser.helpButtonToolTipText=Help
%                     
%       Almost all Swing widgets can be customize this way. You can
%       examine the Swing sources to get these values or check
%       http://www.gargoylesoftware.com/papers/plafdiff.html for
%       a list of them.

