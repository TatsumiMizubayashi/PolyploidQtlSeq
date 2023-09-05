using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Bulk1ディレクトリ
    /// </summary>
    public class Bulk1Directory : ISampleDirectory
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "b1";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "bulk1";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Bulk1 Directory.";

        /// <summary>
        /// Bulk1ディレクトリを作成する。
        /// </summary>
        /// <param name="dirPath">Bulk1ディレクトリのPath</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public Bulk1Directory(string dirPath, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Path = OptionValue.GetValue(LONG_NAME, dirPath, parameterDictionary, userOptionDictionary);
            if (string.IsNullOrWhiteSpace(Path)) throw new ArgumentException(null, nameof(dirPath));

            SampleName = System.IO.Path.GetFileName(Path);
        }

        /// <summary>
        /// Bulk1ディレクトリを作成する。
        /// </summary>
        /// <param name="dirPath">Bulk1ディレクトリのPath</param>
        public Bulk1Directory(string dirPath)
        {
            if (string.IsNullOrWhiteSpace(dirPath)) throw new ArgumentException(null, nameof(dirPath));

            Path = System.IO.Path.GetFullPath(dirPath);
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
