using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.Mapping
{
    /// <summary>
    /// Bulk2ディレクトリ
    /// </summary>
    public class Bulk2Directory : ISampleDirectory
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "b2";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "bulk2";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Bulk2 Directory.";

        /// <summary>
        /// Bulk2ディレクトリを作成する。
        /// </summary>
        /// <param name="dirPath">Bulk2ディレクトリのPath</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public Bulk2Directory(string dirPath, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Path = OptionValue.GetValue(LONG_NAME, dirPath, parameterDictionary, userOptionDictionary);
            if (string.IsNullOrWhiteSpace(Path)) throw new ArgumentException(null, nameof(dirPath));

            SampleName = System.IO.Path.GetFileName(Path);
        }

        public string Path { get; }

        public string SampleName { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Path}";
        }
    }
}
