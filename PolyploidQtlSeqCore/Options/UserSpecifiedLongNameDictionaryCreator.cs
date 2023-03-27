using McMaster.Extensions.CommandLineUtils;

namespace PolyploidQtlSeqCore.Options
{
    /// <summary>
    /// 起動引数情報からユーザーが指定したオプションのLongName辞書を作成する
    /// </summary>
    internal static class UserSpecifiedLongNameDictionaryCreator
    {
        /*
         * 注意：McMaster.Extensions.CommandLineUtilsのVerの違いによって挙動が異なる。
         * Ver3ではユーザー指定の有無はHasValue()で判断できたがVer4ではできない。
         * 
         * Ver4ではオプション指定がない場合はValue()とDefaultValueにコンストラクタで指定した値が入る。
         * ユーザーがオプションを指定した場合はValue()に指定した値が入る。DefaultValueは空になる。
         * 
         */

        /// <summary>
        /// ユーザーが指定したオプションのLongName辞書を作成する。
        /// </summary>
        /// <param name="options">オプションリスト</param>
        /// <returns>ユーザー指定LongName辞書</returns>
        public static IReadOnlyDictionary<string, bool> Create(IReadOnlyCollection<CommandOption> options)
        {
            // LongNameの有無を検索しやすいようにDictionaryで作成する。Valueの型はなんでもいい
            var longNameDictionary = new Dictionary<string, bool>();

            var userSpecifiedOptionQuery = options
                .Where(x => x.DefaultValue != x.Value());

            foreach (var option in userSpecifiedOptionQuery)
            {
                if (string.IsNullOrWhiteSpace(option.LongName)) throw new ArgumentException("LongName is blank.");
                longNameDictionary[option.LongName] = true;
            }

            return longNameDictionary;
        }
    }
}
