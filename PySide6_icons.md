# Icons

## Enum Constant,Visual Description

    QStyle.StandardPixmap.SP_FileIcon,A generic blank document.
    QStyle.StandardPixmap.SP_DirIcon,A standard folder icon.
    QStyle.StandardPixmap.SP_DirOpenIcon,An open folder icon.
    QStyle.StandardPixmap.SP_DirHomeIcon,"A folder with a ""home"" house symbol."
    QStyle.StandardPixmap.SP_DriveHDIcon,A hard drive (what you were looking for earlier).
    QStyle.StandardPixmap.SP_TrashIcon,A wastebasket/delete icon.

## Enum Constant,Visual Description

    QStyle.StandardPixmap.SP_ArrowUp,Upward pointing arrow.
    QStyle.StandardPixmap.SP_ArrowBack,Leftward pointing arrow (Back).
    QStyle.StandardPixmap.SP_ArrowForward,Rightward pointing arrow (Forward).
    QStyle.StandardPixmap.SP_BrowserReload,A circular refresh arrow.
    QStyle.StandardPixmap.SP_TitleBarMaxButton,"The ""Maximize"" window square."

## Enum Constant,Visual Description

    QStyle.StandardPixmap.SP_MessageBoxInformation,"A blue ""i"" circle."
    QStyle.StandardPixmap.SP_MessageBoxWarning,A yellow exclamation triangle.
    QStyle.StandardPixmap.SP_MessageBoxCritical,"A red ""X"" circle."
    QStyle.StandardPixmap.SP_MessageBoxQuestion,A blue question mark.
    QStyle.StandardPixmap.SP_DialogHelpButton,A help/question mark icon.

## Enum Constant,Usage

    QStyle.StandardPixmap.SP_DialogSaveButton,The classic diskette/save icon.
    QStyle.StandardPixmap.SP_DialogOkButton,"A green checkmark or ""OK"" symbol."
    QStyle.StandardPixmap.SP_DialogApplyButton,"Usually a checkmark, used for ""Apply changes."""
    QStyle.StandardPixmap.SP_DialogDiscardButton,"Used for ""Don't Save"" or discarding changes."
    QStyle.StandardPixmap.SP_DialogCancelButton,"A red ""X"" or ""Cancel"" symbol."

## Enum Constant,Visual Representation,Best Used For

    QStyle.StandardPixmap.SP_FileDialogDetailedView,Looks like a list/gear menu,General settings or preferences.
    QStyle.StandardPixmap.SP_ComputerIcon,A monitor/computer,System-wide or hardware settings.
    QStyle.StandardPixmap.SP_CommandLink,An arrow icon,Navigating to a settings sub-page.
    QStyle.StandardPixmap.SP_VistaShield,A security shield,Admin or security settings (Windows).
    QStyle.StandardPixmap.SP_MediaPlay
    QStyle.StandardPixmap.SP_MediaStop

## Constant,Look,Best Use Case

    QStyle.StandardPixmap.SP_FileDialogDetailedView,List with small lines,General preferences or configuration.
    QStyle.StandardPixmap.SP_ComputerIcon,A Monitor,"Hardware, system, or ""Core"" settings."
    QStyle.StandardPixmap.SP_DriveHDIcon,A Hard Drive,Storage or project path settings.
    QStyle.StandardPixmap.SP_TitleBarMenuButton,Small menu/lines,"A ""Hamburger"" style settings menu."
    QStyle.StandardPixmap.SP_FileDialogContentsView,Small grid,Layout or view settings.

## Feature,QStyle Constant

    Save,SP_DialogSaveButton
    Settings,SP_FileDialogDetailedView
    Open/Folder,SP_DirIcon
    Refresh,SP_BrowserReload
    Help,SP_DialogHelpButton
    Home,SP_DirHomeIcon

## Use

btn.setIcon(
    QApplication.style().standardIcon(QStyle.SP_BrowserReload)
).setToolTip("Refresh")
