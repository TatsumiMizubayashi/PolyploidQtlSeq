using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Parent1ディレクトリ
    /// </summary>
    public class Parent1Directory : ISampleDirectory
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "p1";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "parent1";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Parent1 Directory.";

        /// <summary>
        /// Parent1ディレクトリを作成する。
        /// </summary>
        /// <param name="dirPath">Parent1ディレクトリのPath</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public Parent1Directory(string dirPath, IReadOnlyDictionary<string, string> parameterDictionary,
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
