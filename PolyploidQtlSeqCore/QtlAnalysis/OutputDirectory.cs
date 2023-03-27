﻿using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// QTL解析結果出力ディレクトリ
    /// </summary>
    public class OutputDirectory
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "o";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "outputDir";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Output Directory.";

        /// <summary>
        /// 出力ディレクトリを作成する。
        /// </summary>
        /// <param name="outputDirPath">出力ディレクトリPath</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public OutputDirectory(string outputDirPath, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Path = OptionValue.GetValue(LONG_NAME, outputDirPath, parameterDictionary, userOptionDictionary);

            if (string.IsNullOrWhiteSpace(Path)) throw new ArgumentException("Specify the output directory in the -o option.");
        }

        /// <summary>
        /// 出力ディレクトリPATHを取得する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// 出力ディレクトリ内のファイルPathを作成する。
        /// </summary>
        /// <param name="fileNames">ファイル名</param>
        /// <returns>出力ディレクトリ内ファイルPATH</returns>
        internal string CreateFilePath(params string[] fileNames)
        {
            var items = new[] { Path }.Concat(fileNames).ToArray();
            return System.IO.Path.Combine(items);
        }

        /// <summary>
        /// サブディレクトリを作成し、そのPathを返す。
        /// </summary>
        /// <param name="subDirName">サブディレクトリ名</param>
        /// <returns>サブディレクトリのPATH</returns>
        internal string CreateSubDir(string subDirName)
        {
            Create();

            var subDirPath = System.IO.Path.Combine(Path, subDirName);
            if (Directory.Exists(subDirPath)) return subDirPath;

            Directory.CreateDirectory(subDirPath);

            return subDirPath;
        }

        /// <summary>
        /// 出力ディレクトリを作成する。
        /// </summary>
        internal void Create()
        {
            if (Directory.Exists(Path)) return;

            Directory.CreateDirectory(Path);
        }


        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns></returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Path}";
        }
    }
}
