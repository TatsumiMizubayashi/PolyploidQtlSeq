﻿using PolyploidQtlSeqCore.QtlAnalysis.Chr;

namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// BAMファイル
    /// </summary>
    internal class BamFile
    {

        /// <summary>
        /// BAMファイルを作成する。
        /// </summary>
        /// <param name="sampleName">サンプル名</param>
        /// <param name="bamFilePath">BAMファイルPath</param>
        public BamFile(string sampleName, string bamFilePath)
        {
            ArgumentException.ThrowIfNullOrEmpty(sampleName);
            ArgumentException.ThrowIfNullOrWhiteSpace(bamFilePath);

            SampleName = sampleName;
            Path = bamFilePath;
            IndexFilePath = Path + ".bai";
        }

        /// <summary>
        /// サンプル名を取得する。
        /// </summary>
        public string SampleName { get; }

        /// <summary>
        /// BAMファイルのPathを取得する。
        /// </summary>
        public string Path { get; }

        /// <summary>
        /// BAMファイルのIndexFilePathを取得する。
        /// (Index Fileが存在しなくても値を持っている）
        /// </summary>
        private string IndexFilePath { get; }

        /// <summary>
        /// Indexファイルを作成する。
        /// </summary>
        /// <returns></returns>
        public async ValueTask CreateIndexFileAsync()
        {
            await Samtools.CreateIndexFileAsync(Path);
        }

        /// <summary>
        /// Index Fileを持っているかどうかを調査する。
        /// </summary>
        /// <returns>Indexファイルが存在するならtrue</returns>
        public bool ExistsIndexFile()
        {
            return File.Exists(IndexFilePath);
        }

        /// <summary>
        /// BAMファイルとIndexファイルを削除する。
        /// </summary>
        public void Delete()
        {
            if (File.Exists(Path)) File.Delete(Path);
            if (File.Exists(IndexFilePath)) File.Delete(IndexFilePath);
        }

        /// <summary>
        /// BAMヘッダー情報に変換する。
        /// </summary>
        /// <returns>BAMヘッダー</returns>
        public BamHeader ToHeader()
        {
            var t = Samtools.HeaderAsync(this).AsTask();
            t.Wait();

            return t.Result;
        }
    }
}
