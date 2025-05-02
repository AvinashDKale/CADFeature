#include "FileDialog.h"
#include <iostream>


FileDialog::FileDialog() : QFileDialog()
{
	m_btnOpen = NULL;
	m_listView = NULL;
	m_treeView = NULL;
	m_selectedFiles.clear();

	this->setOption(QFileDialog::DontUseNativeDialog, true);
	this->setFileMode(QFileDialog::Directory);
	QList<QPushButton*> btns = this->findChildren<QPushButton*>();
	for (int i = 0; i < btns.size(); ++i) {
		QString text = btns[i]->text();
		if (text.toLower().contains("open") || text.toLower().contains("choose"))
		{
			m_btnOpen = btns[i];
			break;
		}
	}

	if (!m_btnOpen) return;

	m_btnOpen->installEventFilter(this);
	//connect(m_btnOpen, SIGNAL(changed()), this, SLOT(btnChanged()))
	m_btnOpen->disconnect(SIGNAL(clicked()));
	connect(m_btnOpen, SIGNAL(clicked()), this, SLOT(chooseClicked()));


	m_listView = findChild<QListView*>("listView");
	if (m_listView) {
		m_listView->setSelectionMode(QAbstractItemView::ExtendedSelection);
	}

	m_treeView = findChild<QTreeView*>();
	if (m_treeView) {
		m_treeView->setSelectionMode(QAbstractItemView::ExtendedSelection);
	}

	auto a = selectedFiles();
	for (auto i : a) {
		std::cout << "." << i.toStdString() << std::endl;
	}
}

bool FileDialog::eventFilter(QObject* watched, QEvent* event)
{
	QPushButton* btn = qobject_cast<QPushButton*>(watched);
	if (btn)
	{
		if (event->type() == QEvent::EnabledChange) {
			if (!btn->isEnabled()) {
				btn->setEnabled(true);
			}
		}
	}

	return QWidget::eventFilter(watched, event);
}


void FileDialog::chooseClicked()
{
	QModelIndexList indexList = m_listView->selectionModel()->selectedIndexes();
	foreach(QModelIndex index, indexList)
	{
		if (index.column() == 0)
		{
			m_selectedFiles.append(this->directory().absolutePath() + "/" + index.data().toString());
		}
	}
	
	QDialog::accept();
}

QStringList FileDialog::selectedFiles()
{
	return m_selectedFiles;
}
