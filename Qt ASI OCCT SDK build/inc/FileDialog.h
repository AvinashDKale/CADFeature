#pragma once
//Subclass QFileDialog for customize allow select both file/folder

//#include "main.moc"
#include <QWidget>
#include <QCheckBox>
#include <QColorDialog>
#include <QErrorMessage>
#include <QFileDialog>
#include <QFontDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QToolBox>
#include<QListView>
#include<QTreeView>

#include <QGuiApplication>
#include <QStyleHints>
#include <QString>

class FileDialog : public QFileDialog {
	Q_OBJECT
private:
	QListView* m_listView;
	QTreeView* m_treeView;
	QPushButton* m_btnOpen;
	QStringList m_selectedFiles;

public slots:
	void chooseClicked();
public:
	FileDialog();
	QStringList selectedFiles();
	bool eventFilter(QObject* watched, QEvent* event);
};